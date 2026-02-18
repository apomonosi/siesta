# This file implements a strategy to find a PEXSI library compatible
# with the legacy native PEXSI interface in Siesta

# Possible PEXSI providers are ELSI (if compiled with PEXSI)
# and stand-alone PEXSI libraries >= 2.0

   # This might apply only to 2.0, since 2.1 discovery is through
   # library archives and not through a target that might have some
   # properties set

   if (SIESTA_WITH_MPI)
    if (NOT TARGET MPI::MPI_CXX)
      # Find MPI
      find_package(MPI OPTIONAL_COMPONENTS CXX)
      if( NOT MPI_CXX_FOUND )
        message(FATAL_ERROR "MPI_CXX, possibly needed by PEXSI, could not be found by CMake")
      endif()
    endif()
   endif()

   set(_found False)

   if (TARGET elsi::elsi)    # We have requested ELSI support
     if (ELSI_WITH_PEXSI)
       add_library(pexsi-stub INTERFACE)
       target_link_libraries(pexsi-stub INTERFACE elsi::elsi)
       message(STATUS "PEXSI library code for native interface taken from ELSI library")
       add_library(${PROJECT_NAME}::pexsi ALIAS pexsi-stub)
       set(_found True)
     else()
       message(STATUS "ELSI library is available, but it does not offer PEXSI")
     endif()
   else()
      # Maybe search for ELSI,
      # in case it is available and we do not have a stand-alone PEXSI library
   endif()

   # We want to support pre-release PEXSI-2.1, which has some issues
   # with target installation. Hence we provide a "Custom" finder.
   # If we are dealing with PEXSI-2.0, a proper cmake package is
   # available

   if (NOT _found)
     # Fallback to searching for 2.1 package
     message(STATUS "Searching for 2.1 PEXSI library...")
     find_package(CustomPEXSI)
     if (CustomPEXSI_FOUND)
       # This does not seem to be needed, but YMMV
       target_link_libraries(PEXSI::PEXSI INTERFACE MPI::MPI_CXX)
       add_library(${PROJECT_NAME}::pexsi ALIAS PEXSI::PEXSI)
       message(STATUS "PEXSI library taken from a pre-compiled PEXSI>=2.1 installation")
       set(_found True)
     endif()
   endif()

   if (NOT _found)
     # Search for 2.0 cmake package (the only one so far which is sane enough)
     message(STATUS "Searching for 2.0 PEXSI library (cmake config package)...")
     find_package(PEXSI CONFIG QUIET)
     if (PEXSI_FOUND)
       # This does not seem to be needed, but YMMV
       target_link_libraries(PEXSI::PEXSI INTERFACE MPI::MPI_CXX)
       add_library(${PROJECT_NAME}::pexsi ALIAS PEXSI::PEXSI)
       message(STATUS "PEXSI library taken from a pre-compiled PEXSI 2.0 installation")
       set(_found True)
     endif()
   endif()
   

   if (NOT _found)
     message(FATAL_ERROR "PEXSI not found")
   endif()




