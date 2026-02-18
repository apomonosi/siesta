#
#  If you have veclibfort, and it works, you could uncomment these lines
#
set(BLAS_LIBRARY "-lveclibfort" CACHE STRING "blas library chosen")
set(LAPACK_LIBRARY "-lveclibfort" CACHE STRING "lapack library chosen")
#
#  More general settings for use with shell modules
#
# set(LAPACK_LIBRARY "$ENV{LAPACK_LIBS}" CACHE STRING "lapack library chosen")
set(SCALAPACK_LIBRARY "$ENV{SCALAPACK_LIBS}" CACHE STRING "scalapack library chosen")
#
# These can be overridden from the command line
#
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "build_type")
set(Fortran_FLAGS_DEBUG  "-g -O0" CACHE STRING "Fortran debug flags")

#  set(LAPACK_DETECTION OFF CACHE BOOL "LAPACK detection")

