# Credit goes to the DFTB+ developers. The idea is borrowed from their
# amazing CMake recepies.
#
# Toolchain file for
#
# Generic build environment
#
# This is a generic template which probably will not work on your system out of the box. You should
# either modify it or override the variables via command line options to make it to
# work. Alternatively, have a look at the specialized toolchain files in this folder as they may
# give you a better starting point for your build environment.
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
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
  CACHE STRING "Fortran compiler flags for Release with debug info build")

set(Fortran_FLAGS_MINRELSIZE "${CMAKE_Fortran_FLAGS_MINRELSIZE}"
  CACHE STRING "Fortran compiler flags for minimal size build")

set(Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
  CACHE STRING "Fortran compiler flags for Debug build")

set(Fortran_FLAGS_CHECK "${CMAKE_Fortran_FLAGS_CHECK}"
  CACHE STRING "Fortran compiler flags for Debug + checking build")

#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
  CACHE STRING "C compiler flags for Debug build")

