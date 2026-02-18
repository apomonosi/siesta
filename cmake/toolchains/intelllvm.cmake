#
# Toolchain file for
#
# Intel compiler, MKL library. LLVM-based compilers (ifx, icx).
#
# Notes:
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
  set(_host_flags "-xHost")
endif()


#
# Fortran compiler settings
#
set(Fortran_FLAGS_RELEASE "-O2 ${_host_flags} -fp-model=strict"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release with debug info build")

set(Fortran_FLAGS_MINRELSIZE "-Os ${_host_flags} -fp-model=strict"
    CACHE STRING "Fortran compiler flags for minimum size build")

set(Fortran_FLAGS_DEBUG "-g -O0 -traceback"
  CACHE STRING "Fortran compiler flags for Debug build")

set(Fortran_FLAGS_CHECK "-g -O0 -traceback -check all"
  CACHE STRING "Fortran compiler flags for Debug (checking) build")

#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O2 ${_host_flags}"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_RELWITHDEBINFO "-g ${C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for Release with debug info build")

set(C_FLAGS_MINSIZEREL "-Os ${_host_flags}"
  CACHE STRING "C compiler flags for minimum size build")

set(C_FLAGS_DEBUG "-g -O0"
  CACHE STRING "C compiler flags for Debug build")

