#
# Toolchain file for
#
# Craycompiler
#
# Notes:
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually
#

# Specific host variables
# (Currently we don't have them for Cray)
set(_host_flags "")


#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

# O3 is equivalent to fast
set(Fortran_FLAGS_RELEASE "-O2 ${_host_flags}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release with debug info build")

set(Fortran_FLAGS_MINSIZEREL "${_host_flags}"
  CACHE STRING "Fortran compiler flags for minimum size build")

set(Fortran_FLAGS_DEBUG "-g -O0"
  CACHE STRING "Fortran compiler flags for Debug build")

set(Fortran_FLAGS_CHECK "-g -O0"
  CACHE STRING "Fortran compiler flags for Check build")


#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O2 ${_host_flags}"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_RELWITDEBINFO "-g ${C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for RelWithDebInfo build")

set(C_FLAGS_MINSIZEREL "${_host_flags}"
  CACHE STRING "C compiler flags for minimum size build")

set(C_FLAGS_DEBUG "-g -O0"
  CACHE STRING "C compiler flags for Debug build")

