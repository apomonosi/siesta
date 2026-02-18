#
# Vega supercomputer at Maribor, Slovenia, with ELPA GPU support
#
# Using the gompic easyBuild toolchain (GCC, OpenBLAS, OpenMPI, CUDA)
#
# Set up the modules as:
#
# ml gompic/2020b
# ml OpenBLAS/0.3.12-GCC-10.2.0
# ml ScaLAPACK/2.1.0-gompic-2020b
# ml netCDF-Fortran/4.5.3-gompic-2020b
# ml CMake/3.18.4-GCCcore-10.2.0
# ml FFTW/3.3.8-gompic-2020b
#
# OpenBLAS is threaded, so use OpenMP, but set OMP_NUM_THREADS to 1 in most cases.
#
set(SIESTA_WITH_OPENMP "ON" CACHE BOOL "with OpenMP")
#
# To compile ELPA with GPU support, use the following shell script (remove '#-- ')
#
#  Notes: '--enable-nvidia-gpu', and not '--enable-nvidia-sm80-gpu' is used,
#         to avoid kernel naming changes. Support for this needs to be extended.
#
#         The compiler flags are those suggested by the ELPA documentation.
#
#         To run Siesta, make sure that the ELPA lib directory is in the
#         LD_LIBRARY_PATH environment variable.
# -------
#-- LAPACK_ROOT=/cvmfs/sling.si/modules/el7/software/OpenBLAS/0.3.12-GCC-10.2.0
#-- SCALAPACK_ROOT=/cvmfs/sling.si/modules/el7/software/ScaLAPACK/2.1.0-gompic-2020b

#-- FC=mpifort CC=mpicc CXX=mpicxx CPP="gfortran -E" \
#-- FCFLAGS="-O3 -march=native -mavx2 -mfma" \
#-- CFLAGS="-O3 -march=native -mavx2 -mfma  -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize" \
#--  LDFLAGS="-L${SCALAPACK_ROOT}/lib -lscalapack -L${LAPACK_ROOT}/lib -lopenblas -lpthread -lm -ldl" \
#--   ../configure \
#--              --enable-nvidia-gpu --with-NVIDIA-GPU-compute-capability=sm_80 \
#--              --enable-c-tests=no --disable-avx512 --prefix=$HOME/lib/elpa
# -------
#
set(SIESTA_WITH_ELPA "ON" CACHE BOOL "with ELPA")

set(LAPACK_LIBRARY
      "-L /cvmfs/sling.si/modules/el7/software/OpenBLAS/0.3.12-GCC-10.2.0/lib -lopenblas -lpthread -lm -ldl"
        CACHE STRING "lapack library chosen")

set(SCALAPACK_LIBRARY
     "-L/cvmfs/sling.si/modules/el7/software/ScaLAPACK/2.1.0-gompic-2020b/lib -lscalapack"
     CACHE STRING "scalapack library chosen")

set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "build_type")
