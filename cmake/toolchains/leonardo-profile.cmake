##- Toolchain file for use in Leonardo at CINECA, Italy
##-
##- Before you use this file, set the environment modules as follows in the shell.
##- (You can cut the following, put it in a toolchain.sh file, and do
##-      . toolchain.sh
##- )
##-  (These or similar versions for the GNU toolchain)

# module load  gcc/11.3.0
# module load  openblas/0.3.21--gcc--11.3.0
# module load  openmpi/4.1.4--gcc--11.3.0-cuda-11.8
# module load  netlib-scalapack/2.2.0--openmpi--4.1.4--gcc--11.3.0
# module load  netcdf-c/4.9.0--openmpi--4.1.4--gcc--11.3.0
# module load  netcdf-fortran/4.6.0--openmpi--4.1.4--gcc--11.3.0
# module load  libxc/6.2.2--gcc--11.3.0-cuda-11.8
# module load  fftw/3.3.10--openmpi--4.1.4--gcc--11.3.0
# module load  cuda/11.8
##-
#
# You might need to load modules for ELSI, ELPA, etc, as needed
#

set(CMAKE_BUILD_TYPE "Profile" CACHE STRING "build_type")
set(Fortran_FLAGS_PROFILE  "-g -O2 -fno-omit-frame-pointer -funwind-tables" CACHE STRING "Fortran profile flags")

# This seems to be needed because Openblas is mistakenly used in place of Scalapack...

set(SCALAPACK_LIBRARY "-L$ENV{NETLIB_SCALAPACK_HOME} -lscalapack" CACHE STRING "scalapack library")

set(SIESTA_WITH_PROFILE_NVTX "ON" CACHE STRING "use nvtx markers")
set(SIESTA_PROFILE_NVTX_LIBRARY
    "$ENV{CUDA_HOME}/targets/x86_64-linux/lib/libnvToolsExt.so"
    CACHE STRING "nvtx library")



