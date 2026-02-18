##- Toolchain file for use in Marconi100 at CINECA, Italy
##-
##- Before you use this file, set the environment modules as follows in the shell.
##- (You can cut the following, put it in a toolchain.sh file, and do
##-      . toolchain.sh
##- )
##-
##- ml python   # For ninja and meson
##- ml cmake
##- ml gnu/8.4.0
##- ml spectrum_mpi/10.4.0--binary
##- ml blas/3.8.0--gnu--8.4.0
##- ml lapack/3.9.0--gnu--8.4.0
##- ml scalapack/2.1.0--spectrum_mpi--10.4.0--binary      # or maybe 10.3.1 ???
##- ml essl/6.2.1--binary
##- # For the extra bits of LD_LIBRARY_PATH needed for essl... found out with: load("xl/16.1.1--binary")
##- export LD_LIBRARY_PATH="/cineca/prod/opt/compilers/xl/16.1.1/binary/lib:${LD_LIBRARY_PATH}"
##- #
##- ml cuda/11.0
##- export LD_LIBRARY_PATH="/cineca/prod/opt/compilers/cuda/11.0/none/lib64:${LD_LIBRARY_PATH}"
##-

set(LAPACK_LIBRARY
      "-L/cineca/prod/opt/libraries/lapack/3.9.0/gnu--8.4.0/lib -llapack -L/cineca/prod/opt/libraries/essl/6.2.1/binary/lib64 -lessl"
        CACHE STRING "lapack library chosen")
	
set(SCALAPACK_LIBRARY
     "-L/cineca/prod/opt/libraries/scalapack/2.1.0/spectrum_mpi--10.3.1--binary/lib -lscalapack"
     CACHE STRING "scalapack library chosen")

set(SIESTA_PROFILE_NVTX_LIBRARY
    "/cineca/prod/opt/compilers/cuda/11.0/none/targets/ppc64le-linux/lib/libnvToolsExt.so"
    CACHE STRING "nvtx library")
    
set(SIESTA_WITH_NETCDF "OFF" CACHE STRING "netcdf support")    # Do not know yet how to activate it properly
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "build_type")
set(Fortran_FLAGS_RELEASE  "-g -O2" CACHE STRING "Fortran release flags")

#set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "build_type")
#set(Fortran_FLAGS_DEBUG  "-g -O0" CACHE STRING "Fortran debug flags")


