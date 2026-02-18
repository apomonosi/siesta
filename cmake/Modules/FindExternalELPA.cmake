# FindExternalELPA.cmake
#
# Search for ELPA library using pkg-config and extract configuration
# from installed header files.
#
# Input variables (possibly set by caller before find_package):
#   ELPA_PREFER_OPENMP            - If ON, search for elpa_openmp first
#   ELPA_API_VERSION              - User override for ELPA API version (cache)
#   ELPA_GPU_STRING               - User override for GPU string (cache)
#   ELPA_REAL_GPU_KERNEL          - User override for real GPU kernel (cache)
#   ELPA_COMPLEX_GPU_KERNEL       - User override for complex GPU kernel (cache)
#
# Output variables (set by this module):
#   ELPA_FOUND                    - True if ELPA was found
#   ELPA_INCLUDE_DIRS             - Include directories for ELPA
#   ELPA_Fortran_INCLUDE_DIRS     - Fortran module directories
#   ELPA_LIBRARIES                - ELPA libraries
#   ELPA_LINK_LIBRARIES           - Libraries to link against
#   ELPA_HAS_OPENMP               - True if OpenMP version was found
#   ELPA_API_VERSION              - ELPA API version (cache)
#   ELPA_AVAILABLE_REAL_KERNELS   - Available real kernels
#   ELPA_AVAILABLE_COMPLEX_KERNELS - Available complex kernels
#   ELPA_HAS_GPU_SUPPORT          - True if GPU support is available
#   ELPA_GPU_STRING               - GPU string for ELPA 'set' method (cache)
#   ELPA_REAL_GPU_KERNEL          - Real GPU kernel name or -1 (cache)
#   ELPA_COMPLEX_GPU_KERNEL       - Complex GPU kernel name or -1 (cache)
#
# This module also creates an imported target:
#   Elpa::elpa
#
# Non-standard pkg-config issues:
#
# 1. The pkg-config filename might be versioned (e.g., elpa-2020.05.001.rc1.pc)
#    instead of simply 'elpa.pc'. This module searches for common versions.
#
# 2. The include directory entry in the pkg-config file may be missing
#    the 'modules' leaf component. This is handled automatically.

find_package(PkgConfig REQUIRED QUIET)

message(STATUS "Searching for ELPA library")

# If user provided LIB_PATHS, append possible pkgconfig directories
if(DEFINED LIB_PATHS)
  foreach(_lib_path IN LISTS LIB_PATHS)
    if(EXISTS "${_lib_path}/pkgconfig")
      set(ENV{PKG_CONFIG_PATH} "${_lib_path}/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
  endforeach()
endif()

# Initialize OpenMP flag
set(ELPA_HAS_OPENMP FALSE)

# Try to find ELPA via pkg-config
if(PKG_CONFIG_FOUND)
  # If OpenMP is preferred, try elpa_openmp first
  if(ELPA_PREFER_OPENMP)
    message(STATUS "Searching for OpenMP-enabled ELPA (elpa_openmp)")
    pkg_check_modules(ELPA QUIET elpa_openmp)
    
    if(ELPA_FOUND)
      set(ELPA_HAS_OPENMP TRUE)
      message(STATUS "Found elpa_openmp")
    else()
      message(STATUS "elpa_openmp not found, trying standard ELPA")
    endif()
  endif()
  
  # If not found yet (either OpenMP not preferred or elpa_openmp not found)
  if(NOT ELPA_FOUND)
    # Try generic name first
    pkg_check_modules(ELPA QUIET elpa)

    # If not found, try versioned names
    if(NOT ELPA_FOUND)
      pkg_search_module(ELPA QUIET
        elpa-2025.06.001
        elpa-2025.01.002
        elpa-2025.01.001
        elpa-2024.05.001
        elpa-2024.03.001
        elpa-2023.11.001
        elpa-2023.05.001
        elpa-2022.11.001
        elpa-2022.05.001
        elpa-2021.11.002
        elpa-2021.11.001
        elpa-2021.05.002
        elpa-2021.05.001
        elpa-2020.11.001
        elpa-2020.05.001
      )
    endif()
  endif()
endif()

if(NOT ELPA_FOUND)
  message(STATUS "ELPA not found")
  return()
endif()

message(STATUS "ELPA found via pkg-config")
message(STATUS "  ELPA_HAS_OPENMP: ${ELPA_HAS_OPENMP}")
message(STATUS "  ELPA_LIBDIR: ${ELPA_LIBDIR}")
message(STATUS "  ELPA_LIBRARIES: ${ELPA_LIBRARIES}")
message(STATUS "  ELPA_LINK_LIBRARIES: ${ELPA_LINK_LIBRARIES}")
message(STATUS "  ELPA_INCLUDEDIR: ${ELPA_INCLUDEDIR}")
message(STATUS "  ELPA_INCLUDE_DIRS: ${ELPA_INCLUDE_DIRS}")

# Retrieve Fortran compilation flags
pkg_get_variable(ELPA_FCFLAGS ${ELPA_MODULE_NAME} fcflags)

# Fix non-standard setting of Fortran module directory
if(ELPA_FCFLAGS)
  message(STATUS "Extracting Fortran include directories from ELPA FCFLAGS")
  message(STATUS "ELPA_FCFLAGS: ${ELPA_FCFLAGS}")
  set(ELPA_Fortran_INCLUDE_DIRS "")
  
  foreach(flag ${ELPA_FCFLAGS})
    if(flag MATCHES "^-I")
      string(REPLACE "-I" "" cleaned_flag "${flag}")
      list(APPEND ELPA_Fortran_INCLUDE_DIRS "${cleaned_flag}")
    endif()
  endforeach()
else()
  message(STATUS "Adding Fortran include directory manually")
  set(ELPA_Fortran_INCLUDE_DIRS "${ELPA_INCLUDE_DIRS}/modules")
endif()

message(STATUS "  ELPA Fortran include paths: ${ELPA_Fortran_INCLUDE_DIRS}")

# Extract ELPA configuration from installed header files
set(VERSION_CHECK_DIR "${PROJECT_BINARY_DIR}/elpa_version_check_build")
file(MAKE_DIRECTORY ${VERSION_CHECK_DIR})

# Write C source file to extract configuration
file(WRITE "${VERSION_CHECK_DIR}/elpa_config_check.c" [=[
#include <stdio.h>
#include <elpa/elpa_version.h>
#include <elpa/elpa_constants.h>

// Try to include GPU configuration options if available
#ifdef __has_include
#  if __has_include(<elpa/elpa_configured_options.h>)
#    include <elpa/elpa_configured_options.h>
#  else
#    warning "elpa_configured_options.h not found, using defaults (GPU=0)"
#    define ELPA_WITH_NVIDIA_GPU_VERSION 0
#    define ELPA_WITH_AMD_GPU_VERSION 0
#    define ELPA_WITH_SYCL_GPU_VERSION 0
#    define ELPA_MISSING_CONFIG_HEADER 1
#  endif
#else
#  define ELPA_WITH_NVIDIA_GPU_VERSION 0
#  define ELPA_WITH_AMD_GPU_VERSION 0
#  define ELPA_WITH_SYCL_GPU_VERSION 0
#endif

// Ensure macros exist even if header didn't define them
#ifndef ELPA_WITH_NVIDIA_GPU_VERSION
#  define ELPA_WITH_NVIDIA_GPU_VERSION 0
#endif
#ifndef ELPA_WITH_AMD_GPU_VERSION
#  define ELPA_WITH_AMD_GPU_VERSION 0
#endif
#ifndef ELPA_WITH_SYCL_GPU_VERSION
#  define ELPA_WITH_SYCL_GPU_VERSION 0
#endif

int main(void) {
#ifdef ELPA_MISSING_CONFIG_HEADER
    fprintf(stderr, "WARNING: elpa_configured_options.h not found. Using an older ELPA version.\n");
#endif

    printf("API_VERSION=%d\n", ELPA_API_VERSION);
    printf("NVIDIA_GPU=%d\n", ELPA_WITH_NVIDIA_GPU_VERSION);
    printf("AMD_GPU=%d\n", ELPA_WITH_AMD_GPU_VERSION);
    printf("SYCL_GPU=%d\n", ELPA_WITH_SYCL_GPU_VERSION);

    printf("AVAILABLE_REAL_KERNELS=");
    #define X(name, val, available, ...) if (available) printf("%s;", #name);
    ELPA_FOR_ALL_2STAGE_REAL_KERNELS(X)
    printf("\n");
    
    printf("AVAILABLE_COMPLEX_KERNELS=");
    ELPA_FOR_ALL_2STAGE_COMPLEX_KERNELS(X)
    printf("\n");

    return 0;
}
]=])

# Create CMakeLists.txt for the configuration check subproject
file(WRITE "${VERSION_CHECK_DIR}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.6)
project(elpa_config_check C)
add_executable(elpa_config_check elpa_config_check.c)
target_include_directories(elpa_config_check PRIVATE ${ELPA_INCLUDE_DIRS})
]=])

# Configure the version checker
execute_process(
  COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}"
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DELPA_INCLUDE_DIRS=${ELPA_INCLUDE_DIRS}
    ${VERSION_CHECK_DIR}
  WORKING_DIRECTORY ${VERSION_CHECK_DIR}
  RESULT_VARIABLE CONFIG_RESULT
  OUTPUT_QUIET
  ERROR_VARIABLE CONFIG_ERROR
)

if(NOT CONFIG_RESULT EQUAL 0)
  message(WARNING "Failed to configure ELPA config check build: ${CONFIG_ERROR}")
  return()
endif()

# Build the version checker
execute_process(
  COMMAND ${CMAKE_COMMAND} --build .
  WORKING_DIRECTORY ${VERSION_CHECK_DIR}
  RESULT_VARIABLE BUILD_RESULT
  OUTPUT_QUIET
  ERROR_VARIABLE BUILD_ERROR
)

if(NOT BUILD_RESULT EQUAL 0)
  message(WARNING "Failed to build ELPA config check program: ${BUILD_ERROR}")
  return()
endif()

# Run the config check program
execute_process(
  COMMAND ${VERSION_CHECK_DIR}/elpa_config_check${CMAKE_EXECUTABLE_SUFFIX}
  WORKING_DIRECTORY ${VERSION_CHECK_DIR}
  RESULT_VARIABLE RUN_RESULT
  OUTPUT_VARIABLE CONFIG_OUTPUT
  ERROR_VARIABLE CONFIG_STDERR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(RUN_RESULT EQUAL 0)
  # Process output line by line
  # Use a rare character (§) as temporary semicolon replacement
  string(REPLACE ";" "§" TEMP_OUTPUT "${CONFIG_OUTPUT}")
  string(REGEX REPLACE "\n" ";" CONFIG_LINES "${TEMP_OUTPUT}")
  
  foreach(line ${CONFIG_LINES})
    if(line MATCHES "^([^=]+)=(.*)$")
      set(KEY ${CMAKE_MATCH_1})
      string(REPLACE "§" ";" VALUE "${CMAKE_MATCH_2}")

      # Special handling for API_VERSION - allow user override
      if(KEY STREQUAL "API_VERSION")
        if(NOT DEFINED ELPA_API_VERSION)
          set(ELPA_API_VERSION "${VALUE}" CACHE STRING "ELPA API version")
          message(STATUS "  Auto-detected ELPA_API_VERSION = ${VALUE}")
        else()
          message(STATUS "  Using user-defined ELPA_API_VERSION = ${ELPA_API_VERSION}")
        endif()
	
      else()
      
        # All other variables remain regular (non-cache) - these are outputs

        set(ELPA_${KEY} "${VALUE}")
      
        # Display kernel lists nicely
        if(KEY MATCHES "AVAILABLE.*KERNELS")
          message(STATUS "  ELPA ${KEY}:")
          foreach(kernel IN LISTS ELPA_${KEY})
            message(STATUS "    ${kernel}")
          endforeach()
        else()
          message(STATUS "  ELPA ${KEY} = ${VALUE}")
        endif()
	
      endif()
      
    endif()
  endforeach()
  
  # Display any warnings from stderr
  if(CONFIG_STDERR)
    message(STATUS "${CONFIG_STDERR}")
  endif()
else()
  message(WARNING "Failed to run ELPA config check program")
  return()
endif()

# Detect GPU support and automatically set GPU-related variables
set(ELPA_HAS_GPU_SUPPORT FALSE)

# Check if any GPU kernels are available
foreach(kernel IN LISTS ELPA_AVAILABLE_COMPLEX_KERNELS ELPA_AVAILABLE_REAL_KERNELS)
  string(FIND "${kernel}" "GPU" _has_gpu_kernel)
  if(NOT _has_gpu_kernel EQUAL -1)
    set(ELPA_HAS_GPU_SUPPORT TRUE)
    break()
  endif()
endforeach()

# Function to set GPU kernel cache variables based on kernel name
# Only sets if not already defined by user
function(set_gpu_kernel kernel)
  string(FIND "${kernel}" "_REAL_" _is_real)
  string(FIND "${kernel}" "_COMPLEX_" _is_complex)
  
  if(NOT _is_real EQUAL -1)
    if(NOT DEFINED ELPA_REAL_GPU_KERNEL)
      set(ELPA_REAL_GPU_KERNEL "${kernel}" CACHE STRING "ELPA real GPU kernel")
      message(STATUS "  Auto-detected ELPA_REAL_GPU_KERNEL: ${kernel}")
    endif()
  elseif(NOT _is_complex EQUAL -1)
    if(NOT DEFINED ELPA_COMPLEX_GPU_KERNEL)
      set(ELPA_COMPLEX_GPU_KERNEL "${kernel}" CACHE STRING "ELPA complex GPU kernel")
      message(STATUS "  Auto-detected ELPA_COMPLEX_GPU_KERNEL: ${kernel}")
    endif()
  endif()
endfunction()

if(ELPA_HAS_GPU_SUPPORT)
  message(STATUS "ELPA has GPU support")
  
  # Report any user-defined values
  if(DEFINED ELPA_REAL_GPU_KERNEL)
    message(STATUS "  User-defined ELPA_REAL_GPU_KERNEL: ${ELPA_REAL_GPU_KERNEL}")
  endif()
  if(DEFINED ELPA_COMPLEX_GPU_KERNEL)
    message(STATUS "  User-defined ELPA_COMPLEX_GPU_KERNEL: ${ELPA_COMPLEX_GPU_KERNEL}")
  endif()
  if(DEFINED ELPA_GPU_STRING)
    message(STATUS "  User-defined ELPA_GPU_STRING: ${ELPA_GPU_STRING}")
  endif()
  
  # Auto-detect GPU vendor and set appropriate kernel variables
  foreach(kernel IN LISTS ELPA_AVAILABLE_REAL_KERNELS ELPA_AVAILABLE_COMPLEX_KERNELS)
    # Stop if everything is already defined
    if(DEFINED ELPA_GPU_STRING AND 
       DEFINED ELPA_REAL_GPU_KERNEL AND 
       DEFINED ELPA_COMPLEX_GPU_KERNEL)
      break()
    endif()
    
    string(FIND "${kernel}" "_NVIDIA_GPU" _is_nvidia)
    string(FIND "${kernel}" "_NVIDIA_SM80_GPU" _is_nvidia_sm80)
    string(FIND "${kernel}" "_AMD_GPU" _is_amd)
    string(FIND "${kernel}" "_INTEL_GPU_SYCL" _is_intel)
    string(FIND "${kernel}" "_COMPLEX_GPU" _is_generic)
    
    if(NOT _is_nvidia_sm80 EQUAL -1)
      set_gpu_kernel("${kernel}")
      if(NOT DEFINED ELPA_GPU_STRING)
        set(ELPA_GPU_STRING "nvidia-gpu" CACHE STRING "ELPA GPU string for 'set' method")
        message(STATUS "  Auto-detected ELPA_GPU_STRING: nvidia-gpu (SM80)")
      endif()
    elseif(NOT _is_nvidia EQUAL -1)
      set_gpu_kernel("${kernel}")
      if(NOT DEFINED ELPA_GPU_STRING)
        set(ELPA_GPU_STRING "nvidia-gpu" CACHE STRING "ELPA GPU string for 'set' method")
        message(STATUS "  Auto-detected ELPA_GPU_STRING: nvidia-gpu")
      endif()
    elseif(NOT _is_amd EQUAL -1)
      set_gpu_kernel("${kernel}")
      if(NOT DEFINED ELPA_GPU_STRING)
        set(ELPA_GPU_STRING "amd-gpu" CACHE STRING "ELPA GPU string for 'set' method")
        message(STATUS "  Auto-detected ELPA_GPU_STRING: amd-gpu")
      endif()
    elseif(NOT _is_intel EQUAL -1)
      set_gpu_kernel("${kernel}")
      if(NOT DEFINED ELPA_GPU_STRING)
        set(ELPA_GPU_STRING "intel-gpu" CACHE STRING "ELPA GPU string for 'set' method")
        message(STATUS "  Auto-detected ELPA_GPU_STRING: intel-gpu")
      endif()
    elseif(NOT _is_generic EQUAL -1)
      set_gpu_kernel("${kernel}")
      if(NOT DEFINED ELPA_GPU_STRING)
        set(ELPA_GPU_STRING "gpu" CACHE STRING "ELPA GPU string for 'set' method")
        message(STATUS "  Auto-detected ELPA_GPU_STRING: gpu (generic)")
      endif()
    endif()
  endforeach()
else()
  message(STATUS "ELPA does not have GPU support")
  if(NOT DEFINED ELPA_GPU_STRING)
    set(ELPA_GPU_STRING "no-gpu" CACHE STRING "ELPA GPU string for 'set' method")
  endif()
  if(NOT DEFINED ELPA_REAL_GPU_KERNEL)
    set(ELPA_REAL_GPU_KERNEL -1 CACHE STRING "No ELPA real GPU kernel")
  endif()
  if(NOT DEFINED ELPA_COMPLEX_GPU_KERNEL)
    set(ELPA_COMPLEX_GPU_KERNEL -1 CACHE STRING "No ELPA complex GPU kernel")
  endif()
endif()

# Create imported target
add_library(Elpa::elpa INTERFACE IMPORTED)
target_link_libraries(Elpa::elpa INTERFACE ${ELPA_LINK_LIBRARIES})
target_include_directories(Elpa::elpa INTERFACE 
  ${ELPA_Fortran_INCLUDE_DIRS}
  ${ELPA_INCLUDE_DIRS}
)

message(STATUS "ELPA configuration complete")
