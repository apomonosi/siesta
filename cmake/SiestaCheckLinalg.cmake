# Run sanity checks on the libraries.
# In particular, check that the return convention is appropriate on MacOS
# (See https://github.com/mcg1969/vecLibFort for more details)
#

include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

set(CMAKE_REQUIRED_QUIET OFF)


message(STATUS "Checking that BLAS library works...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Check that sgemm can be found (blas check)
if( TARGET BLAS::BLAS )
  set(CMAKE_REQUIRED_LIBRARIES BLAS::BLAS)
elseif( TARGET LAPACK::LAPACK )
  set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
endif()

#
# Check linking 
#
check_fortran_source_compiles(
"
external :: sgemm
call sgemm()
end
"
  blas_has_sgemm SRC_EXT F90)

if (NOT blas_has_sgemm)
  message(WARNING
    "---------------------------------------------\n"
    " BLAS library cannot link properly\n"
    " Please check the library linking string found or used by CMake\n"
    "---------------------------------------------")
  message(FATAL_ERROR " BLAS library does not link properly")
endif()

#
# Check that BLAS/LAPACK actually works
#
siesta_try_run_default(blas_sgemm_basic TRUE)
check_fortran_source_runs(
"
real :: a(4), b(4), c(4)
external :: sgemm
! Initialize 2x2 matrices (stored column-wise)
! A = [[1,2],[3,4]], B = identity matrix
a = [1.0, 3.0, 2.0, 4.0]
b = [1.0, 0.0, 0.0, 1.0] 
c = [0.0, 0.0, 0.0, 0.0]
! Perform C = A * B
call sgemm('N', 'N', 2, 2, 2, 1.0, a, 2, b, 2, 0.0, c, 2)
end
"
blas_sgemm_basic SRC_EXT F90)

if (NOT blas_sgemm_basic)
  message(WARNING
    "---------------------------------------------\n"
    " Cannot run simple program using the BLAS library's sgemm \n"
    " Please check the library linking string found or used by CMake\n"
    "---------------------------------------------")
  message(FATAL_ERROR " Cannot run simple program using BLAS library")
endif()

#
# Now check the return-value convention
#
siesta_try_run_default(blas_cdotu_return_convention TRUE)
check_fortran_source_runs(
"
complex :: c
complex, dimension(2) :: a = [ 1, 2 ]
complex :: cdotu
external :: cdotu
c=cdotu(2,a(:),1,a(:),1)
end
"
blas_cdotu_return_convention SRC_EXT F90)
unset(CMAKE_REQUIRED_LIBRARIES)

if (NOT blas_cdotu_return_convention)
  message(WARNING
    "---------------------------------------------\n"
    " BLAS library uses wrong return-value convention!"
    " This is likely to happen on MacOS if the default Accelerate framework is used"
    " You can install veclibfort (https://github.com/mcg1969/vecLibFort)"
    " and then set either of:\n"
    "   -DBLAS_LIBRARY=-lveclibfort"
    " or"
    "   -DLAPACK_LIBRARY=-lveclibfort\n"
    " in your cmake invocation."
    " Where the key point is to have the BLAS library linked first."
    " Alternatively you can install OpenBLAS/BLIS and set the variables accordingly.\n"
    "---------------------------------------------")
  message(FATAL_ERROR " BLAS library uses wrong return-value convention!!!")
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)


### CHECK LAPACK

message(STATUS "Checking that LAPACK library works...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Test for lapack
if(TARGET LAPACK::LAPACK)
  set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
endif()
check_fortran_source_compiles(
"
external :: dsysv
call dsysv()
end
"
lapack_has_dsysv SRC_EXT F90)

if (NOT lapack_has_dsysv)
  message(WARNING
    "---------------------------------------------\n"
    " LAPACK library cannot link properly"
    " Please check the library linking string found or used by CMake\n"
    "---------------------------------------------")
  message(FATAL_ERROR "  *** LAPACK library does not link properly")
endif()


unset(CMAKE_REQUIRED_LIBRARIES)

list(POP_BACK CMAKE_MESSAGE_INDENT)


### CHECK ScaLAPACK
if(NOT SIESTA_WITH_MPI)
  return()
endif()

message(STATUS "Checking that ScaLAPACK library works...")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Test for scalapack
if(TARGET SCALAPACK::SCALAPACK)
  set(CMAKE_REQUIRED_LIBRARIES SCALAPACK::SCALAPACK)
endif()
check_fortran_source_compiles(
"
external :: blacs_gridinit
call blacs_gridinit()
end
"
scalapack_has_blacs_gridinit SRC_EXT F90)

if (NOT scalapack_has_blacs_gridinit)
  message(WARNING
    "---------------------------------------------\n"
    " ScaLAPACK library cannot link properly"
    " Please check the library linking string found or used by CMake"
    "---------------------------------------------")
  message(FATAL_ERROR "  *** ScaLAPACK library does not link properly")
endif()
unset(CMAKE_REQUIRED_LIBRARIES)

list(POP_BACK CMAKE_MESSAGE_INDENT)

