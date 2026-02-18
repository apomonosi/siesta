#
# Toolchain file for
#
# GNU compiler
#
# Notes:
#
#  * Settings here should work out of the box on Ubuntu (tested on 18.4). Other build environments
#    may need some fine tuning.
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually
#

# Specific host variables
if( CMAKE_CROSSCOMPILING OR NOT SIESTA_WITH_HOST_OPTIMIZATION )
  set(_host_flags "")
else()
  set(_host_flags "-march=native")
endif()

#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O3 ${_host_flags}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE} -fbacktrace"
  CACHE STRING "Fortran compiler flags for Release with debug info build")

set(Fortran_FLAGS_MINSIZEREL "-Os ${_host_flags}"
  CACHE STRING "Fortran compiler flags for minimum size build")

set(Fortran_FLAGS_DEBUG "-g -Og -fbacktrace"
  CACHE STRING "Fortran compiler flags for Debug build")

set(Fortran_FLAGS_CHECK "-g -Og -fbacktrace -fcheck=all"
  CACHE STRING "Fortran compiler flags for Debug build")


#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O3 ${_host_flags}"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_RELWITDEBINFO "-g ${C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for RelWithDebInfo build")

set(C_FLAGS_MINSIZEREL "-Os ${_host_flags}"
  CACHE STRING "C compiler flags for minimum size build")

set(C_FLAGS_DEBUG "-g -Wall -pedantic"
  CACHE STRING "C compiler flags for Debug build")

