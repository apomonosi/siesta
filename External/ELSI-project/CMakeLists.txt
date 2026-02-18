
list(APPEND CMAKE_MESSAGE_CONTEXT elsi)

if(ELSI_ALLOW_PRECOMPILED)

  message(STATUS " ** Searching for pre-compiled ELSI library...")
  message(STATUS " ** ... You might need to define extra variables.")
  message(STATUS " ** ... Continue only if you know what you are doing.")

  # We need to know:
  #  -- whether ELSI uses an external ELPA library
  #  -- whether ELSI has GPU support

  list(APPEND CMAKE_MESSAGE_CONTEXT search-for-elsi)
  include(search_for_elsi)
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)

  if (ELSI_FOUND)
    message(STATUS " ** ... found. Note that features available might vary")

    if (ELSI_WITH_PEXSI AND (NOT ELSI_HAS_PEXSI))
      message(STATUS "** Warning: Installed ELSI library does not have PEXSI")
      message(STATUS "** Warning: Falling back to compiling it from sub-project.")
      set(ELSI_FOUND FALSE)
    else()
      list(POP_BACK CMAKE_MESSAGE_CONTEXT)
      return()
    endif()
  endif()

endif()

# -- Use ELSI as a sub-project

# C and C++ standards applying
# Older standards tend to be more forgiving with some
# constructs. This is needed in the scotch library, for example,
# with modern versions of gcc.
# This is a more portable way to specify the
# standard level (instead of, for example -std=c99)
# 
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# To override internal option() statements with variables
# set by the parent project
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
#
# To accept <PACKAGE>_ROOT (uppercase version) as hint
#
set(CMAKE_POLICY_DEFAULT_CMP0144 NEW)

# Note that ELPA2024 support has been removed from ELSI
# The 2021 version has been retained, but its use in Siesta has not been tested much,
# and it does not work on the Mac.
# It might not work for GPU either (missing .cu files -- to be reported upstream)
option(ELSI_WITH_ELPA2021 "Use 2021 ELPA in ELSI" FALSE)
if(APPLE AND ELSI_WITH_ELPA2021)
   message(FATAL_ERROR "ELPA2021 does not work on the Mac. Remove ELSI_WITH_ELPA2021.")
endif()

# Should these options be passed to the top level?
# Note that ELSI_WITH_PEXSI triggers the enabling of CXX, which can only be
# done at the top level. It is set to FALSE here.
#
option(ELSI_WITH_GPU "Use GPU in ELSI" FALSE)
option(ELSI_WITH_EXTERNAL_ELPA "Use an external ELPA library in ELSI" FALSE)
option(ELSI_WITH_PEXSI "Enable PEXSI in ELSI" FALSE)
option(ELSI_WITH_EXTERNAL_PEXSI "Use an external PEXSI library in ELSI" FALSE)
option(ELSI_WITH_TESTS "Compile (Fortran) tests in ELSI" TRUE)
option(ELSI_WITH_C_TESTS "Compile (C) tests in ELSI" FALSE)
option(ELSI_INSTALL_PKGCONFIG_FILE "Install pkgconfig file" FALSE)

# The ELSI subproject will also honor the variables for compiler options:
#
#  ELSI_Fortran_FLAGS
#  ELSI_C_FLAGS
#  ELSI_CXX_FLAGS
#
#  which will be inserted in the CMAKE_Lang_FLAGS_CONFIG variables, where
#  CONFIG is the (uppercase) value of CMAKE_BUILD_CONFIG.
#  (or the custom ELSI_CUSTOM if the build configuration is not defined)

# For GPU (currently only CUDA), the variables honored are:
#
# CMAKE_CUDA_ARCHITECTURES   (defaults to 80 (A100 class))
# CMAKE_CUDA_FLAGS           (see default below)
#

# Setting variables for ELSI subproject:

# So that we do not need to use the MPI compilers
# If this feature ever stops working, we can enable the following

if (FALSE)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} CACHE STRING "new Fortran compiler")
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE STRING "new C compiler")
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "new CXX compiler")
else()
  set(ELSI_USE_FIND_PACKAGE_MPI ON)
endif()

# Setup proper dependency targets list for ELSI
set(ELSI_EXTERNAL_TARGETS)

# Note that when using ELSI_USE_FIND_PACKAGE_MPI, ELSI's CMake scripts add transitive
# dependencies to MPI::MPI_C and MPI::MPI_Fortran

list(APPEND ELSI_EXTERNAL_TARGETS SCALAPACK::SCALAPACK LAPACK::LAPACK)

#
# Options for built-in libs
#

# 2020 and 2021 are the only supported versions in ELSI as of now
# Actually, 2020 is recommended
if(ELSI_WITH_ELPA2021)
  message(WARNING "ELPA_2021 is not recommended. Falling back to 2020")
  ###  set(USE_ELPA_2021 ON)
endif()


if(ELSI_WITH_PEXSI)
  set(ENABLE_PEXSI ON)
  set(ELSI_HAS_PEXSI ON)
endif()
#
# Direct transcription of wrapper options to ELSI options via variables
#
set(ENABLE_TESTS ${ELSI_WITH_TESTS})
set(ENABLE_C_TESTS ${ELSI_WITH_C_TESTS})


# Handle external solver libraries in ELSI (just ELPA for now)

if(ELSI_WITH_EXTERNAL_ELPA)
  if (NOT TARGET Elpa::elpa)
    message(STATUS "Looking for pre-compiled ELPA library for ELSI sub-project...")
    ## If this ELSI wrapper is ever made portable to other projects, we might need to deal with OpenMP here
    ## As of now, the ELPA_PREFER_OPENMP option is set in the parent, and this block to search for
    ## the ELPA target is actually not needed
    find_package(ExternalELPA)
    if (NOT ELPA_FOUND)
      message(WARNING "External ELPA library not found for ELSI sub-project. Falling back to internal ELPA in ELSI")
    endif()
  endif()

  if (TARGET Elpa::elpa)
    message(STATUS "Using external ELPA for ELSI compilation")
    list(APPEND ELSI_EXTERNAL_TARGETS Elpa::elpa)
    set(USE_EXTERNAL_ELPA ON)
    set(ELSI_HAS_EXTERNAL_ELPA ON)
  else()
    set(ELSI_HAS_EXTERNAL_ELPA OFF)
    set(ELPA_API_VERSION "20200000")
  endif()
endif()

# Logic for GPU support

if(ELSI_HAS_EXTERNAL_ELPA)
  if(ELPA_HAS_GPU_SUPPORT)
    message(STATUS "External ELPA library in use has GPU support (${ELPA_GPU_STRING})")
    set(ELSI_ELPA_HAS_GPU_SUPPORT TRUE)
    set(ELSI_ELPA_GPU_STRING ${ELPA_GPU_STRING})
    set(ELSI_ELPA_REAL_GPU_KERNEL ${ELPA_REAL_GPU_KERNEL})
    set(ELSI_ELPA_COMPLEX_GPU_KERNEL ${ELPA_COMPLEX_GPU_KERNEL})
  else()
    message(STATUS "External ELPA library in use has NO GPU support")
    # We could also error out here
    if(ELSI_WITH_GPU)
      message(WARNING "External ELPA library in use for ELSI has NO GPU support")
    endif()
    set(ELSI_ELPA_HAS_GPU_SUPPORT FALSE)
    set(ELSI_ELPA_GPU_STRING "no-gpu")
    set(ELSI_ELPA_REAL_GPU_KERNEL -1)
    set(ELSI_ELPA_COMPLEX_GPU_KERNEL -1)
  endif()
else(ELSI_HAS_EXTERNAL_ELPA)
  if(ELSI_WITH_GPU)
    # We will compile ELSI with GPU support. Currently, only CUDA
    set(USE_GPU_CUDA ON)

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
    message(STATUS "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")

    if(NOT DEFINED CMAKE_CUDA_FLAGS)
      set(CMAKE_CUDA_FLAGS "-arch=sm_${CMAKE_CUDA_ARCHITECTURES} -O3 -std=c++11")
    endif()
    message(STATUS "CMAKE_CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")

    find_package(CUDAToolkit REQUIRED)
    list(APPEND ELSI_EXTERNAL_TARGETS CUDA::cublasLt CUDA::cublas CUDA::cudart)
    if(ELSI_WITH_ELPA2021)
      set(ELSI_ELPA_GPU_STRING "nvidia-gpu")
      message(STATUS "** Using ELPA2021 within ELSI with GPU support")
      message(STATUS "** Note that it might not have all the needed files...")
    else()
      set(ELSI_ELPA_GPU_STRING "gpu")
      message(STATUS "Using ELPA2020 within ELSI with GPU support")
    endif()

    set(ELSI_HAS_CUDA ON)
    set(ELSI_ELPA_HAS_GPU_SUPPORT TRUE)
    set(ELSI_ELPA_REAL_GPU_KERNEL "ELPA_2STAGE_REAL_GPU")
    set(ELSI_ELPA_COMPLEX_GPU_KERNEL "ELPA_2STAGE_COMPLEX_GPU")
    message(STATUS "ELSI_ELPA_HAS_GPU_SUPPORT: ${ELSI_ELPA_HAS_GPU_SUPPORT}")
  else(ELSI_WITH_GPU)
    set(ELSI_HAS_CUDA OFF)
    set(ELSI_ELPA_HAS_GPU_SUPPORT FALSE)
    set(ELSI_ELPA_GPU_STRING "no-gpu")
    set(ELSI_ELPA_REAL_GPU_KERNEL -1)
    set(ELSI_ELPA_COMPLEX_GPU_KERNEL -1)
  endif(ELSI_WITH_GPU)
endif(ELSI_HAS_EXTERNAL_ELPA)

message(STATUS " ** External targets passed to ELSI:")

foreach(_tgt ${ELSI_EXTERNAL_TARGETS})
  message(STATUS " --> ${_tgt}")
endforeach()

# Pass this to the LIBS variable, which in ELSI will actually end up holding target names
# Note that it is important that this variable is not a cache variable in ELSI

set(LIBS ${ELSI_EXTERNAL_TARGETS})

include(FetchContent)

if(EXISTS "${PROJECT_SOURCE_DIR}/External/ELSI-project/elsi_interface/CMakeLists.txt")
  set(FETCHCONTENT_SOURCE_DIR_ELSI "${PROJECT_SOURCE_DIR}/External/ELSI-project/elsi_interface")
  message(STATUS "... Configuring existing ELSI sources in ${PROJECT_SOURCE_DIR}/External/ELSI-project/elsi_interface")
else()
  message(STATUS "... Fetching ELSI sources from Git repository")
endif()

FetchContent_Declare(
  elsi
  GIT_REPOSITORY https://gitlab.com/elsi_project/elsi_interface.git
  # Keep this in sync with submodule commit reference
  GIT_TAG v2.12.0
  ## GIT_TAG  master
  ## For later -- OVERRIDE_FIND_PACKAGE
  )

FetchContent_MakeAvailable(elsi)

set(ELSI_FOUND TRUE)

if (NOT TARGET elsi::elsi)
  message(STATUS "-- Creating elsi::elsi target")
  add_library(elsi::elsi ALIAS elsi)
endif()

# Pass result variables upstream

set(ELSI_FOUND ${ELSI_FOUND} PARENT_SCOPE)
set(ELSI_HAS_PEXSI ${ELSI_HAS_PEXSI} PARENT_SCOPE)
set(ELSI_HAS_EXTERNAL_ELPA ${ELSI_HAS_EXTERNAL_ELPA} PARENT_SCOPE)
set(ELSI_HAS_CUDA ${ELSI_HAS_CUDA} PARENT_SCOPE)
set(ELSI_ELPA_HAS_GPU_SUPPORT ${ELSI_ELPA_HAS_GPU_SUPPORT} PARENT_SCOPE)
set(ELPA_HAS_GPU_SUPPORT ${ELPA_HAS_GPU_SUPPORT} PARENT_SCOPE)
set(ELSI_ELPA_GPU_STRING ${ELSI_ELPA_GPU_STRING} PARENT_SCOPE)
set(ELSI_ELPA_REAL_GPU_KERNEL ${ELSI_ELPA_REAL_GPU_KERNEL} PARENT_SCOPE)
set(ELSI_ELPA_COMPLEX_GPU_KERNEL ${ELSI_ELPA_COMPLEX_GPU_KERNEL} PARENT_SCOPE)

list(POP_BACK CMAKE_MESSAGE_CONTEXT)

