#
# Leonardo at CINECA with nvfortran
#
# Use the following stack (or similar versions)
# ml nvhpc/23.1
# ml openmpi/4.1.4--nvhpc--23.1-cuda-11.8
# ml openblas/0.3.21--nvhpc--23.1
# ml netlib-scalapack/2.2.0--openmpi--4.1.4--nvhpc--23.1
# ml fftw/3.3.10--openmpi--4.1.4--nvhpc--23.1
# ml netcdf-fortran/4.6.0--openmpi--4.1.4--nvhpc--23.1
# ml cuda/11.8   (Maybe only at runtime, to avoid warnings at compilation)
#
# This is needed so that nvfortran does not interpret "\" in strings as
# an escape character (mostly for DFTD3).
# This could also be added generally if the compiler is detected to be NVHPC
#
# These are now trapped in the main configuration 
##set(Fortran_FLAGS "-Mbackslash" CACHE STRING "Fortran debug flags")
#
# Some compiler issues remain here
#
##set(SIESTA_WITH_FLOOK OFF CACHE BOOL "do not use flook yet -- issues remain")
#
# MPI interfaces have to be deactivated
#
##set(SIESTA_WITH_NO_MPI_INTERFACES ON CACHE BOOL "remove legacy MPI interfaces")

# This seems to be needed because Openblas might be mistakenly used in place of Scalapack...

set(SCALAPACK_LIBRARY "-L$ENV{NETLIB_SCALAPACK_HOME} -lscalapack" CACHE STRING "scalapack library")