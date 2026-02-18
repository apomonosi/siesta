#
#  This file will enable linking to the Accelerate framework in MacOS, which
#  provides BLAS and LAPACK.
#
#  HOWEVER, it is only for test purposes. The return-value convention of the Accelerate framework
#  is that of f2c, not gfortran's. 
#  See https://github.com/mcg1969/vecLibFort for an alternative
#
#  The standard Siesta linear-algebra library checks should catch this.
#
set(SIESTA_WITH_MPI OFF CACHE BOOL "Do not need MPI for this test")
#
if(APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -framework Accelerate")
    # If you encounter symbol issues:
    # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -framework Accelerate -fno-underscoring")
endif()
#


