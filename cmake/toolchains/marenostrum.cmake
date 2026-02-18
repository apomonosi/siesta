#
# Marenostrum with Intel 2017.4 toolchain. NO OPENMP
#
set(SIESTA_WITH_OPENMP "OFF" CACHE BOOL "OpenMP setting")

set(LAPACK_LIBRARY
    "-L/apps/INTEL/2017.4/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
        CACHE STRING "lapack library chosen")
	
set(SCALAPACK_LIBRARY
    "-L/apps/INTEL/2017.4/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
     CACHE STRING "scalapack library chosen")
     

