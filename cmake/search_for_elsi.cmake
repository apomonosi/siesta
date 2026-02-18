#
#
# -- Search for ELSI using pkg-config data by default
#

include(CheckFortranSourceRuns)
include(CheckFortranSourceCompiles)


 # This seems to work better with external libraries (including CUDA), since somehow all the pieces
 # for the link are specified in the elsi.pc file.
 # (see https://gitlab.com/elsi_project/elsi_interface/-/issues/57)

 find_package(PkgConfig REQUIRED)

 pkg_check_modules(ELSI elsi>=2.9.1)

 if(ELSI_FOUND)
   message(STATUS "Found ELSI library (pkg-conf file) Version: ${ELSI_VERSION}")

   message(DEBUG "ELSI Libdir: ${ELSI_LIBDIR}")
   message(DEBUG "ELSI Libraries: ${ELSI_LIBRARIES}")
   message(STATUS "ELSI_LINK_LIBRARIES: ${ELSI_LINK_LIBRARIES}")
   message(STATUS "ELSI_INCLUDEDIR: ${ELSI_INCLUDEDIR}")

   add_library(elsi::elsi INTERFACE IMPORTED GLOBAL)
   target_link_libraries(elsi::elsi INTERFACE  ${ELSI_LINK_LIBRARIES})
   target_include_directories(elsi::elsi INTERFACE ${ELSI_INCLUDEDIR})

 else()

  set(ELSI_MIN_VERSION "2.9.0")
  find_package(elsi ${ELSI_MIN_VERSION} QUIET)
  if (elsi_FOUND)
   message(STATUS "Found ELSI library (CMake package) Version: ${elsi_VERSION}")

   if(NOT TARGET elsi::elpa)
     # ELSI uses an external ELPA library
     # To link it properly, we need to add the right ELPA directory,
     # as the CMake-based ELSI package just links "-lelpa".

     message(STATUS "The found ELSI library needs an external ELPA library")
     message(STATUS "If the following does not work, re-enable SEARCH_ELSI_PKGCONF...")

     if(NOT ELPA_LIBDIR)
       # We have not searched for ELPA...
       # Since Siesta can use a native interface to ELPA, we just suggest to
       # the user to set -DSIESTA_WITH_ELPA=ON. Alternatively, the user should provide
       # the right value in ELPA_LIBDIR
       message(STATUS "**")
       message(STATUS "The found ELSI library needs an external ELPA library")
       message(STATUS "Either set -DSIESTA_WITH_ELPA=ON")
       message(STATUS "  (recommended, and you get also the native Siesta interface to ELPA)")
       message(STATUS "or set -DELPA_LIBDIR=/path/to/elpa/libdir")
       message(STATUS "  (this might not work with all compilers)")
       message(FATAL_ERROR "**")
     else()
       message(STATUS "The found ELSI library needs an external ELPA library")
       message(STATUS "Adding ELPA directory for ELSI linking: ${ELPA_LIBDIR}")
       target_link_directories(elsi::elsi INTERFACE "${ELPA_LIBDIR}")
       #
       # Note: Some compilers (e.g., maybe the NEC one) migh need access to
       # the full chain of modules. In that case we would also need:
       #     target_include_directories(elsi::elsi INTERFACE ${ELPA_FORTRAN_INC_DIRS})
       # It might be better to just 'include(search_for_elpa)' here
       # and use the variables. 
     endif()
   endif()
  endif()
  set(ELSI_FOUND ${elsi_FOUND})
 endif()

 if(NOT ELSI_FOUND)
   return()
 endif()

  set(CMAKE_REQUIRED_LIBRARIES elsi::elsi MPI::MPI_Fortran)

  # Check first whether we can link a very simple program with an "external" call
  # without adding the C++ runtime, just MPI.
  
  check_fortran_source_compiles("program p; call elsi_get_pexsi_enabled(i); end"
    ELSI_LINKS_OK SRC_EXT F90)
  unset(CMAKE_REQUIRED_LIBRARIES)
    
  if (ELSI_LINKS_OK)
    set(CMAKE_REQUIRED_LIBRARIES elsi::elsi MPI::MPI_Fortran)
      check_fortran_source_runs("program p; call elsi_get_pexsi_enabled(i); if (i == 0) ERROR STOP 1; end"
      ELSI_HAS_PEXSI_ALT SRC_EXT F90)
    check_fortran_source_runs("program p; call elsi_get_magma_enabled(i); if (i == 0) ERROR STOP 1; end"
      ELSI_HAS_MAGMA SRC_EXT F90)
    check_fortran_source_runs("program p; call elsi_get_eigenexa_enabled(i); if (i == 0) ERROR STOP 1; end"
      ELSI_HAS_EIGENEXA SRC_EXT F90)
    check_fortran_source_runs("program p; call elsi_get_sips_enabled(i); if (i == 0) ERROR STOP 1; end"
      ELSI_HAS_SIPS SRC_EXT F90)
    unset(CMAKE_REQUIRED_LIBRARIES)

    message(STATUS "Compiled check for PEXSI: ${ELSI_HAS_PEXSI_ALT}")
    set(ELSI_HAS_PEXSI ${ELSI_HAS_PEXSI_ALT})
  else()
     message(WARNING "ELSI library cannot be interrogated easily about features")
     # Must use command-line variables to set
     if(NOT DEFINED CACHE{ELSI_WITH_PEXSI})
       set(ELSI_HAS_PEXSI FALSE)
     endif()
     # These to be implemented if needed (can use the calls in the client program)
     set(ELSI_HAS_MAGMA FALSE)
     set(ELSI_HAS_EIGENEXA FALSE)
     set(ELSI_HAS_SIPS FALSE)
  endif()

if(ELSI_HAS_PEXSI)

  # PEXSI is written in C++, so we need this
  # ... and maybe also to search for MPI_CXX??

  enable_language(CXX)
  # Avoid to link everything associated with PEXSI with C++
  set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
  # We still need to link with the C++ runtime
  # A more portable solution used now involves adding a dummy C++ library to the top link step

endif()

