cmake_minimum_required(VERSION 3.20 FATAL_ERROR)


# CMAKE_MODULE_PATH is locally searched for *.cmake files
# before CMAKE_PREFIX_PATH.
# We will prefer to use CMAKE_MODULE_PATH
list(PREPEND
  CMAKE_MODULE_PATH
  # add the generic configuration cmake directory
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  # add the external modules directory
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules"
  )

# Define project details
# Add [siesta] to the output
list(APPEND CMAKE_MESSAGE_CONTEXT siesta)
project(
  SIESTA
  LANGUAGES Fortran C
  HOMEPAGE_URL "https://siesta-project.org/siesta"
  DESCRIPTION "A first-principles materials simulation code using DFT."
  )

set(SIESTA_AUTHOR  "Siesta group")
set(SIESTA_LICENSE "GPLv3")
#-----

# Use new rpath behavior
cmake_policy(SET CMP0042 NEW)

# Follow GNU conventions for installing directories
# This will use the bin/ include/ etc. folder structure
include(GNUInstallDirs)

# For dependent options
include(CMakeDependentOption)

# Include Siesta specific utility functions
include(SiestaUtils)
include(SiestaConda)
include(SiestaPrintBuildInfo)
siesta_util_ensure_out_of_source_build()

if(PROJECT_IS_TOP_LEVEL)
  if(SIESTA_ALLOW_LEGACY_OPTION_NAMES) # Please try to avoid this
    include(SiestaEnableLegacyOptions)
  else()
    include(SiestaBlockLegacyOptions)
  endif()
else()
  message(STATUS "*** Siesta is being configured as a subproject***")
  message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
  message(STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
  include(SiestaBlockLegacyOptions)
endif()

# Compiler flags
include(SiestaFlags)
include(SiestaCheckFlags)


# -- Options ----------------------------------------

# This is just to have the user choose how to do the build
# generally they don't need to be touched, but it can be useful
# when siesta is used as a submodule.
option(SIESTA_TESTS "Siesta: Build test-suite" ${PROJECT_IS_TOP_LEVEL})
if( SIESTA_TESTS )
  enable_testing()
endif()
option(BUILD_SHARED_LIBS "Siesta: Build using shared libraries" FALSE)
option(SIESTA_SHARED_LIBS "Siesta: Build as a shared library" "${BUILD_SHARED_LIBS}")
if( SIESTA_SHARED_LIBS )
  include(CheckPIESupported)
  check_pie_supported()
  message(STATUS "Enabling position independent code")
  # Ensure BUILD_SHARED_LIBS is set
  set(BUILD_SHARED_LIBS ON CACHE BOOL "Building shared libs" FORCE)
elseif( BUILD_SHARED_LIBS )
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Building shared libs" FORCE)
endif()

option(SIESTA_SET_RPATH "Siesta: add rpath information to executable and libraries" ON)
if( SIESTA_SET_RPATH )
  message(STATUS "Enabling setting of rpath in installed objects")

  # Common settings
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR};${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  else()
    # Linux-specific settings
    set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR};${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    # Use RUNPATH instead of RPATH
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-new-dtags")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--enable-new-dtags")
  endif()

endif()

if(SIESTA_HAS_IPO)
  option(SIESTA_WITH_IPO "Siesta: enable IPO/LTO optimization for Siesta targets" FALSE)
  if( SIESTA_WITH_IPO )
    message(STATUS "Enabling IPO/LTO for targets")
  endif()
elseif(SIESTA_WITH_IPO)
  message(FATAL_ERROR [==[
Siesta detected something wrong IPO for your compilers.

Either your default flags are erroneous or IPO simply does not
work with your compiler...

Please reconfigure without -DSIESTA_WITH_IPO
]==])
endif()


option(SIESTA_INSTALL "Siesta: Install project" ${PROJECT_IS_TOP_LEVEL})

# this should be advanced, as it is primarily useful for users of Siesta as a
# sub-project.
mark_as_advanced(BUILD_SHARED_LIBS)
mark_as_advanced(SIESTA_INSTALL SIESTA_TESTS SIESTA_SHARED_LIBS)

option(SIESTA_WITH_GRID_SP "Siesta: Use single-precision for grid magnitudes" FALSE)

option(SIESTA_BUILD_DOCS "Siesta: Create documentation (WIP)" FALSE)
option(SIESTA_WITH_PEXSI "Siesta: Use the PEXSI library (direct interface)" FALSE)
option(SIESTA_WITH_WANNIER90 "Siesta: Use the wannier90 interface" FALSE)
option(SIESTA_WITH_CHESS "Siesta: Use the CheSS linear-scaling library from BigDFT" FALSE)

# If SIESTA_WITH_UNIT_CONVENTION is set as a normal variable by a
# super-project, the normal variable will take precedence in
# principle, but versions of CMake older than 3.21 will *remove* the
# normal variable when executing the set(...CACHE...) statement...  So
# this must be wrapped

if(PROJECT_IS_TOP_LEVEL)
 # Command-line options take precedence over this, since the 'force' option is not
 # given, for obvious reasons
 set(SIESTA_WITH_UNIT_CONVENTION "CODATA2018"
   CACHE STRING
   "Decide which units ${PROJECT_NAME} should be compiled with [original,legacy,codata2018]"
   )
 # Define the available options for cmake-gui
 set_property(CACHE SIESTA_WITH_UNIT_CONVENTION
   PROPERTY STRINGS
   "original" "legacy" "CODATA2018"
 )
# Do not make this public
 mark_as_advanced(SIESTA_WITH_UNIT_CONVENTION)

else()
 if (NOT DEFINED SIESTA_WITH_UNIT_CONVENTION)
  set(SIESTA_WITH_UNIT_CONVENTION "CODATA2018")
 endif()
endif()

# ------------------------------------------------


# Use Pkg-Config if available
# Many of our sub-inclusions requires this one, so we will catch errors here
# To not confuse users!
find_package(PkgConfig REQUIRED)


# -- Version information
message(CHECK_START "Checking Siesta version")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Before searching for any versions, we create the required version files
# that can be used to track version information
execute_process(COMMAND sh "${PROJECT_SOURCE_DIR}/Tools/version_generate.sh"
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  OUTPUT_VARIABLE SIESTA_VERSION
  RESULT_VARIABLE VERSION_RESULT
  ERROR_VARIABLE VERSION_STDERR)

# fix omission from vgen which sometimes has a newline
string(REGEX REPLACE "\n" " " SIESTA_VERSION "${SIESTA_VERSION}")
string(STRIP "${SIESTA_VERSION}" SIESTA_VERSION)
message(STATUS "SIESTA_VERSION: ${SIESTA_VERSION} (err=${VERSION_RESULT})")
# Pretty print the error content
if(VERSION_RESULT EQUAL 2)
  set(tmp CHECK_FAIL)
  set(tmp_msg "not found")
elseif(VERSION_RESULT EQUAL 3)
  set(tmp CHECK_PASS)
  set(tmp_msg "found development version")
else()
  set(tmp CHECK_PASS)
  set(tmp_msg "found release version")
endif()

message(STATUS "${VERSION_STDERR}")
list(POP_BACK CMAKE_MESSAGE_INDENT)
message(${tmp} "${tmp_msg}")


# Look for blas, if needed, it is not required if LAPACK ships BLAS
find_package(CustomBlas)

# Look for LAPACK
find_package(CustomLapack REQUIRED)

if(NOT TARGET MPI::MPI_Fortran)
 # Find MPI
 find_package(MPI COMPONENTS Fortran)
endif()


# AG: Scalapack-related logic should perhaps be updated

if(MPI_Fortran_FOUND)
  if( NOT DEFINED SIESTA_WITH_MPI)
    find_package(CustomScalapack)
    if( NOT SCALAPACK_FOUND )
      message(WARNING
        "MPI is found, but ScaLAPACK library cannot be found (or compiled against).\n"
        "If parallel support is required please supply the ScaLAPACK library with appropriate "
        "flags:\n"
        " -DSCALAPACK_LIBRARY=<lib>")
      set(MPI_Fortran_FOUND FALSE)
    endif()
  elseif( SIESTA_WITH_MPI )
    find_package(CustomScalapack)
  endif()
endif()

option(SIESTA_WITH_MPI "Build Siesta with MPI support" ${MPI_Fortran_FOUND})

if( SIESTA_WITH_MPI )
  if(NOT MPI_Fortran_FOUND)
    message(WARNING
      "MPI could not be found by CMake"
      "Try setting MPI_Fortran_COMPILER to an MPI wrapper, or the MPI_HOME env variable"
      "If unsuccessful, set -DSIESTA_WITH_MPI=OFF")
    message(FATAL_ERROR "MPI was requested, but could not be found by CMake")
  endif()
  if( NOT SCALAPACK_FOUND )
    message(WARNING
      "ScaLAPACK could not be found by CMake"
      "Try adding SCALAPACK_LIBRARIES with appropriate flags"
      "If ScaLAPACK is not present MPI needs to be turned off:\n"
      "  -DSIESTA_WITH_MPI=OFF"
      )
    message(FATAL_ERROR "MPI was requested, but ScaLAPACK could not be found by CMake")
  endif()


else()
  message(STATUS "MPI is disabled")

  # Add the fakeMPI directory to the include path to ensure it's just working

endif()




# Check for ELSI support
option(SIESTA_WITH_ELSI "Siesta: Use the ELSI library interface" ${SIESTA_WITH_MPI})

# Now we can check for linear algebra stuff
include(SiestaCheckLinalg)
include(SiestaCheckLinalgFeatures)

#---------------------------- OpenMP
# OpenMP *must* be user-requested since it brings additional
# overhead. We should also, for consistency check whether the LAPACK
# and/or BLAS libraries support OpenMP in some form.
# There is little gain when using OpenMP and BLAS libraries are not threaded.
option(SIESTA_WITH_OPENMP "Build with OpenMP support (very few utilities will benefit)" FALSE)
if( SIESTA_WITH_OPENMP )

  if(NOT TARGET OpenMP::OpenMP_Fortran)
    find_package(OpenMP COMPONENTS Fortran C)

    if( OpenMP_FOUND )

      if(NOT OpenMP_Fortran_HAVE_OMPLIB_MODULE)
        message(ERROR_FATAL "OpenMP does not have the 'omp_lib' module used by Siesta OpenMP")
      endif()

    else()
      message(ERROR_FATAL "OpenMP support requested but could not find any OpenMP flags")
    endif()
  endif()

endif()

# NETCDF -- First, a fault check for case-sensitive arguments
# Instead of relying on multiple variants, we should just
# kill the build.
# The problem is that NetCDF *may* require additional
# variables in correct case format. And if we allow upper-case
# for the _ROOT variable, we should do so for all of the variables
# which turns out to be a mess... (we can do a loop if we have everything,
# but for now lets not do this.
foreach(name IN ITEMS "PATH" "ROOT")
  if( (DEFINED NETCDF_${name}) AND (NOT DEFINED NetCDF_${name}) )
    message(WARNING
      "Searching for NetCDF requires the correct case of the search."
      "The build script detected that NETCDF_${name} is defined whereas "
      "NetCDF_${name} is not defined.\n"
      "The search for the NetCDF library requires that you use the "
      "NetCDF_* format for variables. E.g. do this:\n"
      " -DNetCDF_${name}=${NETCDF_${name}}\n"
      "on the command-line instead.")
    message(FATAL_ERROR "Wrong case-notation of NetCDF search variables.")
  endif()
endforeach()

find_package(NetCDF 4.0.0 COMPONENTS Fortran C)

message(DEBUG "NetCDF_Fortran_INCLUDE_DIRS: ${NetCDF_Fortran_INCLUDE_DIRS}")
message(DEBUG "NetCDF_C_INCLUDE_DIRS: ${NetCDF_C_INCLUDE_DIRS}")

if (DEFINED SIESTA_WITH_NetCDF AND (NOT DEFINED SIESTA_WITH_NETCDF))
  message(FATAL_ERROR "Please use SIESTA_WITH_NETCDF and not SIESTA_WITH_NetCDF as the option name")
endif()

option(SIESTA_WITH_NETCDF "Build Siesta with NetCDF support" ${NetCDF_Fortran_FOUND})

if (SIESTA_WITH_NETCDF AND (NOT NetCDF_Fortran_FOUND))
  message(WARNING
    "NetCDF_Fortran is requested but could not be found by CMake.\n"
    "Try setting variables these variables for searching:\n"
    " - NetCDF_ROOT or NetCDF_PATH\n"
    " - NetCDF_INCLUDE_DIR or NetCDF_Fortran_INCLUDE_DIR\n"
    " - NetCDF_INCLUDE_DIRS or NetCDF_Fortran_INCLUDE_DIRS\n"
    "to appropriate variables for discovering NetCDF fortran libraries."
    "Another possibility would be to add the NetCDF cmake package path to\n"
    "  CMAKE_PREFIX_PATH\n"
    "for automatic discovery."
    "If still unsuccessful, please disable NetCDF support with:\n"
    "  -DSIESTA_WITH_NETCDF=OFF\n")
  message(FATAL_ERROR "NetCDF could not be found by CMake")
endif()

option(SIESTA_WITH_NCDF "Use the NCDF library" ${SIESTA_WITH_NETCDF})


if( SIESTA_WITH_NETCDF )

  if (SIESTA_WITH_MPI)
    # Note that MPI must be enabled for this to make sense
    option(SIESTA_WITH_NCDF_PARALLEL "Use NCDF_PARALLEL" ${NetCDF_PARALLEL})
  else()
    option(SIESTA_WITH_NCDF_PARALLEL "Use NCDF_PARALLEL" FALSE)
  endif()

else()

  # We "shadow" any cache values of SIESTA_WITH_NCDF and
  # SIESTA_WITH_NCDF_PARALLEL with the values of normal variables.
  # Avoid modifying any cache values, which might not have the desired
  # effect in all versions of CMake (see policy 0126)

  set(SIESTA_WITH_NCDF FALSE CACHE BOOL "" FORCE)
  set(SIESTA_WITH_NCDF_PARALLEL FALSE CACHE BOOL "" FORCE)

endif()

if(SIESTA_WITH_MPI)

  set(SIESTA_WITH_MPI_INTERFACES "auto" CACHE STRING  "Interface choice undefined")

  # Define the available options
  set_property(CACHE SIESTA_WITH_MPI_INTERFACES
    PROPERTY STRINGS
   "f08" "legacy" "none" "auto"
   )

  if(DEFINED  SIESTA_WITH_MPI_INTERFACES)
    # Add validation
    if(NOT SIESTA_WITH_MPI_INTERFACES MATCHES "^(f08|legacy|none|auto)$")
      message(FATAL_ERROR "SIESTA_WITH_MPI_INTERFACES must be one of: f08, legacy, none, auto (got '${SIESTA_WITH_MPI_INTERFACES}')")
    endif()
  endif()

  if(MPI_Fortran_HAVE_F08_MODULE)
    message(STATUS "MPI Fortran 2008 module (mpi_f08) is available.")
    if (SIESTA_WITH_MPI_INTERFACES STREQUAL "auto")
      message(STATUS "Activating use of mpi_f08 module")
      set(SIESTA_WITH_MPI_INTERFACES "f08" CACHE STRING  "Use mpi_f08 module" FORCE)
    endif()
  else()
    message(STATUS "MPI Fortran 2008 module (mpi_f08) is NOT available.")
    if (SIESTA_WITH_MPI_INTERFACES STREQUAL "f08")
      message(WARNING "Falling back to legacy MPI interfaces")
    endif()
    set(SIESTA_WITH_MPI_INTERFACES "legacy" CACHE STRING  "Use legacy MPI interfaces" FORCE)
  endif()

  # This message might alert the user if a compiler
  # cannot process the legacy interfaces and does not have the f08 ones
  if (SIESTA_WITH_MPI_INTERFACES STREQUAL "legacy")
    message(STATUS "Attempting to use legacy MPI interfaces")
    message(STATUS "Set SIESTA_WITH_MPI_INTERFACES=none if not adequate")
  endif()

endif()

# Handle external dependencies
# Use proper order for xmlf90 and libpsml dependencies
if (NOT TARGET libfdf::libfdf)
  find_package(Customlibfdf REQUIRED)
endif()

if (NOT TARGET xmlf90::xmlf90)
  find_package(Customxmlf90 REQUIRED)
endif()

if (NOT TARGET libpsml::libpsml)
  find_package(Customlibpsml REQUIRED)
endif()
# ensure libpsml links with xmlf90
target_link_libraries(libpsml::libpsml
  INTERFACE
  xmlf90::xmlf90)

if(NOT LIBPSML_USES_PROCEDURE_POINTER )

  siesta_suffix_install()

  # Work around bug in shared library linking
  # This should only temporarily be a hack since future
  # psml versions will always rely on procedure pointers
  siesta_add_library(
    ${PROJECT_NAME}.psml_wrappers
    NO_LINKER_FLAGS
    STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/Src/psml_wrappers.F90"
  )

  siesta_suffix_uninstall()

  # One cannot use OBJECT libraries as their objects are not
  # passed down to subsequent dependent usages. I.e. that would cause
  # name-collisions for static libraries.
  target_link_libraries(libpsml::libpsml
    INTERFACE
    ${PROJECT_NAME}.psml_wrappers
  )

endif()

if (NOT TARGET Libxc::xc_Fortran)
  find_package(CustomLibxc)
  if (TARGET Libxc::xc_Fortran)
    message(VERBOSE "Defining new target libxc::XC_Fortran...")
    add_library(libxc::XC_Fortran INTERFACE IMPORTED)
    target_link_libraries(libxc::XC_Fortran
      INTERFACE
      Libxc::xc_Fortran)
  endif()
endif()

option(SIESTA_WITH_LIBXC "Include libxc support" ${LIBXC_Fortran_FOUND})

if( SIESTA_WITH_LIBXC AND (NOT LIBXC_Fortran_FOUND) )
  message(WARNING
    "Libxc support has been requested ${SIESTA_WITH_LIBXC} : ${LIBXC_Fortran_FOUND}, but the fortran libxc library cannot be found.")
  message(FATAL_ERROR "Libxc could not be found by CMake")
endif()

if (NOT TARGET libgridxc::libgridxc)
  find_package(CustomLibGridxc REQUIRED)
endif()
if(NOT LIBGRIDXC_USES_PROCEDURE_POINTER)

  siesta_suffix_install()

  siesta_add_library(
    ${PROJECT_NAME}.gridxc_wrappers
    NO_LINKER_FLAGS
    STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/Src/gridxc_wrappers.F90"
  )

  siesta_suffix_uninstall()

  target_link_libraries(libgridxc::libgridxc
    INTERFACE
    ${PROJECT_NAME}.gridxc_wrappers
  )

endif()

option(SIESTA_WITH_DFTD3 "Support for DFT-D3 corrections" TRUE)
if(SIESTA_WITH_DFTD3)
  add_subdirectory("External/DFTD3")
  if(NOT TARGET s-dftd3::s-dftd3)
    message(WARNING "Cannot configure DFT-D3 module. Check sources in External/DFTD3")
    #  Shadow the cache variable with a normal variable
    set(SIESTA_WITH_DFTD3 FALSE CACHE BOOL "" FORCE)
  endif()
endif()

# Search for flook ------
# We make it ON by default, except if using the NVHPC compiler
# ... this is a bit convoluted
#
if(CMAKE_Fortran_COMPILER_ID MATCHES NVHPC)
  set(_enable_flook_default False)
  message(STATUS "SIESTA_WITH_FLOOK will be set to OFF by default as some NVHPC compiler versions have trouble with it.")
  message(STATUS "You can still override the setting if it actually works, or if you have a working compiled library.")
else()
  set(_enable_flook_default True)
endif()

siesta_add_subdirectory_option(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/External/Lua-Engine/"
  HELP "Include support for Lua engine through flook (custom MD + SCF interaction)"
  NITEMS 2
  DEFAULT ${_enable_flook_default}
  OPTION SIESTA_WITH_FLOOK
  )
#------

# Search for wannier90 code
if( SIESTA_WITH_WANNIER90 )
  add_subdirectory("External/Wannier")
endif()


if( SIESTA_WITH_MPI )
  # no need to search for ELPA in non-mpi environments

  cmake_dependent_option(ELPA_PREFER_OPENMP
                       "Prefer OpenMP version of ELPA"
                       ON                       # Default when SIESTA_WITH_OPENMP is ON
                       "SIESTA_WITH_OPENMP"     # Condition to check
                       OFF)                     # Forced value when SIESTA_WITH_OPENMP is OFF

  message(STATUS "About to search for ELPA: ELPA_PREFER_OPENMP:  ${ELPA_PREFER_OPENMP}")
  find_package(ExternalELPA)
  option(SIESTA_WITH_ELPA "Include ELPA support" ${ELPA_FOUND})

  if (ELPA_FOUND)
    set(SIESTA_ELPA_HAS_GPU_SUPPORT ${ELPA_HAS_GPU_SUPPORT})
    message(STATUS "SIESTA_ELPA_HAS_GPU_SUPPORT:  ${SIESTA_ELPA_HAS_GPU_SUPPORT}")
    set(SIESTA_ELPA_GPU_STRING ${ELPA_GPU_STRING})
    set(SIESTA_ELPA_REAL_GPU_KERNEL ${ELPA_REAL_GPU_KERNEL})
    set(SIESTA_ELPA_COMPLEX_GPU_KERNEL ${ELPA_COMPLEX_GPU_KERNEL})
    if (ELPA_PREFER_OPENMP AND (NOT ELPA_HAS_OPENMP))
      message(WARNING "Could not find an OpenMP-enabled external ELPA library")
    endif()
  endif()
endif()

# For ELSI-PEXSI, we want it ON by default if ELSI is ON, but need to track if user explicitly set it
if(NOT DEFINED SIESTA_WITH_ELSI_PEXSI)
    set(SIESTA_WITH_ELSI_PEXSI ${SIESTA_WITH_ELSI} CACHE BOOL "Siesta: Compile PEXSI support in ELSI")
    set(ELSI_PEXSI_USER_SET FALSE CACHE INTERNAL "Tracks if PEXSI was explicitly set by user")
else()
    set(ELSI_PEXSI_USER_SET TRUE CACHE INTERNAL "Tracks if PEXSI was explicitly set by user")
endif()

if( SIESTA_WITH_ELSI )

 if(SIESTA_WITH_ELPA)
   set(ELSI_WITH_EXTERNAL_ELPA TRUE)
 endif()

 if(SIESTA_WITH_ELSI_PEXSI)
    find_package(BISON)
    find_package(FLEX 2.6.4)
    if(NOT BISON_FOUND OR NOT FLEX_FOUND OR NOT "${FLEX_VERSION}" STREQUAL "2.6.4")
       if(ELSI_PEXSI_USER_SET)
           # User explicitly wanted PEXSI, so fail if dependencies not found
           if(NOT BISON_FOUND)
              message(FATAL_ERROR "ELSI-PEXSI support explicitly requested but bison not found")
           endif()
           if(NOT FLEX_FOUND)
              message(FATAL_ERROR "ELSI-PEXSI support explicitly requested but flex not found")
           endif()
	   if(NOT "${FLEX_VERSION}" STREQUAL "2.6.4")
              message(FATAL_ERROR "To compile ELSI-PEXSI you need flex 2.6.4")
           endif()
       else()
           # Dependencies not found but user didn't explicitly request PEXSI,
           # so disable it with a warning
           set(SIESTA_WITH_ELSI_PEXSI FALSE CACHE BOOL "Siesta: Compile PEXSI support in ELSI" FORCE)
           message(WARNING "Disabling PEXSI support in ELSI due to missing dependencies:")
           if(NOT BISON_FOUND)
              message(WARNING "  - bison not found")
           endif()
           if(NOT FLEX_FOUND)
              message(WARNING "  - flex not found")
           endif()
	   if(NOT "${FLEX_VERSION}" STREQUAL "2.6.4")
              message(WARNING "  - flex version 2.6.4 not found")
           endif()
       endif()
    endif()
 endif()

 if(SIESTA_WITH_ELSI_PEXSI)
   set(ELSI_WITH_PEXSI TRUE)
   # Languages can only be declared at the top level
   enable_language(CXX)
   set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
 endif()

 # Use subordinate directory context to encapsulate ELSI variables
 # ELSI_FOUND and ELSI_HAS_PEXSI et al are passed upwards from the ELSI context

 list(APPEND CMAKE_MESSAGE_CONTEXT elsi-subproject)
 add_subdirectory("External/ELSI-project")

 if (NOT ELSI_FOUND)
  message(FATAL_ERROR
    "ELSI support has been requested, but the ELSI library cannot be found nor compiled.")
 endif()

 if( ELSI_HAS_PEXSI )
   if (NOT TARGET cxx_dummy_lib)
     # To link properly the C++ runtime, create a small dummy C++ library
     message(STATUS "Compiling dummy C++ library")
     add_subdirectory(Src/stubs)
   endif()
 endif()
 list(POP_BACK CMAKE_MESSAGE_CONTEXT)

 if (SIESTA_WITH_ELPA)
   set(SIESTA_HAS_ELPA_THROUGH_ELSI FALSE)
 else()
   set(SIESTA_HAS_ELPA_THROUGH_ELSI TRUE)
   set(SIESTA_WITH_ELPA TRUE)
   message(STATUS "ELSI can provide ELPA: ${SIESTA_WITH_ELPA}")

   set(SIESTA_ELPA_HAS_GPU_SUPPORT ${ELSI_ELPA_HAS_GPU_SUPPORT})
   message(STATUS "SIESTA_ELPA_HAS_GPU_SUPPORT:  ${SIESTA_ELPA_HAS_GPU_SUPPORT}")
   set(SIESTA_ELPA_GPU_STRING ${ELSI_ELPA_GPU_STRING})
   set(SIESTA_ELPA_REAL_GPU_KERNEL ${ELSI_ELPA_REAL_GPU_KERNEL})
   set(SIESTA_ELPA_COMPLEX_GPU_KERNEL ${ELSI_ELPA_COMPLEX_GPU_KERNEL})
 endif()

endif()

set(SIESTA_NEEDS_ELPA_TARGET FALSE)
set(SIESTA_NEEDS_ELSI_TARGET_FOR_ELPA FALSE)

if(SIESTA_WITH_ELPA)
    if(SIESTA_HAS_ELPA_THROUGH_ELSI)
        set(SIESTA_NEEDS_ELSI_TARGET_FOR_ELPA TRUE)
    else()
        set(SIESTA_NEEDS_ELPA_TARGET TRUE)
    endif()
else()
    set(SIESTA_ELPA_HAS_GPU_SUPPORT FALSE)
    message(STATUS "SIESTA_ELPA_HAS_GPU_SUPPORT:  ${SIESTA_ELPA_HAS_GPU_SUPPORT}")
    set(SIESTA_ELPA_GPU_STRING "no-gpu")
    set(SIESTA_ELPA_REAL_GPU_KERNEL -1)
    set(SIESTA_ELPA_COMPLEX_GPU_KERNEL -1)
endif()

message(STATUS "SIESTA_NEEDS_ELPA_TARGET: ${SIESTA_NEEDS_ELPA_TARGET}")
message(STATUS "SIESTA_NEEDS_ELSI_TARGET_FOR_ELPA: ${SIESTA_NEEDS_ELSI_TARGET_FOR_ELPA}")

# Improve messages here
if( SIESTA_WITH_CHESS )
  find_package(CustomCHESS)
endif()

if( SIESTA_WITH_PEXSI )

   if( NOT SIESTA_WITH_MPI )
     message(FATAL_ERROR "PEXSI support is only available with MPI")
   endif()

   enable_language(CXX)
   #
   # Avoid to link everything associated with PEXSI with C++
   # (i.e., the linker will be the Fortran compiler, and
   # the runtime C++ library will be added to the link statement)
   #
   set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)

   include(search_for_native_PEXSI_provider)

endif(SIESTA_WITH_PEXSI)

# For Poisson and STM/ol-stm.
find_package(FFTW COMPONENTS DOUBLE_LIB)
option(SIESTA_WITH_FFTW "Build with FFTW support." ${FFTW_DOUBLE_LIB_FOUND})


# install the suffix mechanism
siesta_suffix_install()

# The MPI library handles other things, and then returns ASAP!
# The .mpi library is not created.
add_subdirectory("Src/MPI")

add_subdirectory("Src/libsys")
add_subdirectory("Src/libunits")
add_subdirectory("Src/easy-fdict")
add_subdirectory("Src/libxc-compat")
add_subdirectory("Src/ncps/src")
add_subdirectory("Src/psoplib/src")
add_subdirectory("Src/MatrixSwitch/src")
if( SIESTA_WITH_NCDF )
  add_subdirectory("Src/easy-ncdf")
endif()

if( SIESTA_WITH_PROFILE_NVTX )
  add_library(${PROJECT_NAME}.nvtx-c INTERFACE IMPORTED)
  target_link_libraries(${PROJECT_NAME}.nvtx-c INTERFACE ${CMAKE_DL_LIBS})
  set_target_properties(${PROJECT_NAME}.nvtx-c
       PROPERTIES
       INTERFACE_LINK_LIBRARIES
       "${SIESTA_PROFILE_NVTX_LIBRARY}")
endif()

# Siesta sources
add_subdirectory("Src")

# -- Utilities
add_subdirectory("Util")
add_subdirectory("Pseudo")

# Package license files
if( SIESTA_INSTALL )
  install(
    FILES
    "COPYING"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
  )
endif()

# Docs
# Enable SIESTA_BUILD_DOCS to create Doxygen+ford-based documentation (wip)
add_subdirectory(Docs)

# add a few very basic tests for now


if( SIESTA_TESTS )
  include(SiestaTestMPIConfig)

  if( SIESTA_WITH_MPI )
    add_subdirectory("Util/MPI_test")
  endif()

  add_subdirectory(Tests)

endif()


# Condense all build-information
siesta_print_build_info()


# Print information about siesta and libsiesta
siesta_print_start_section("Siesta build")

siesta_print_feature_info(
  NAME ${PROJECT_NAME}
  TARGETS
    ${PROJECT_NAME}::siesta
    ${PROJECT_NAME}::libsiesta
    ${PROJECT_NAME}::libhsx
  VARIABLES
    SIESTA_WITH_UNIT_CONVENTION
    SIESTA_WITH_MPI
    SIESTA_WITH_OPENMP
    SIESTA_WITH_NETCDF
    SIESTA_WITH_NCDF
    SIESTA_WITH_LIBXC
    SIESTA_WITH_PEXSI
    SIESTA_WITH_FFTW
    SIESTA_WITH_DFTD3
    SIESTA_WITH_WANNIER90
    SIESTA_WITH_CHESS
    SIESTA_WITH_ELSI
    SIESTA_WITH_ELPA
    SIESTA_WITH_FLOOK
    # Compilation specific details
    SIESTA_WITH_HOST_OPTIMIZATION
    SIESTA_WITH_IPO
    SIESTA_WITH_NETCDF_PARALLEL
    SIESTA_WITH_NO_MPI_INTERFACES
    SIESTA_INSTALL
    BUILD_SHARED_LIBS
    SIESTA_SHARED_LIBS
    SIESTA_LINKER_FLAGS_PRE
    SIESTA_LINKER_FLAGS_POST
    SIESTA_LINKER_FLAGS
    SIESTA_TESTS
    SIESTA_TESTS_MPI_MIN_NUMPROCS
    SIESTA_TESTS_MPI_NUMPROCS
    SIESTA_TESTS_MPI_MAX_NUMPROCS
    SIESTA_TOOLCHAIN
    SIESTA_SUFFIX
    SIESTA_EXECUTABLE_SUFFIX
    SIESTA_SHARED_LIBRARY_SUFFIX
    SIESTA_STATIC_LIBRARY_SUFFIX
    SIESTA_INTEGER_KINDS # for cross-compilation
    SIESTA_REAL_KINDS # for cross-compilation
    SIESTA_WITH_PROFILE_NVTX
    SIESTA_PROFILE_NVTX_LIBRARY
    SIESTA_BUILD_DOCS
  MSGON
    "The ${PROJECT_NAME} build information will be listed here."
    "In particular how the targets and the variables affecting binaries are interpreted."
    "It can be used to debug how the final linker lines for siesta "
    "is handled."
  REQUIRED
  )

siesta_print_end_section()

siesta_suffix_uninstall()

# Ensure we remove the Siesta context
list(POP_BACK CMAKE_MESSAGE_CONTEXT)
