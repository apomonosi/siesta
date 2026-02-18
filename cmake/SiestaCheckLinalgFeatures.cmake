include(CheckFortranSourceCompiles)

message(STATUS "Checking linear algebra compatibility and extensions")
list(APPEND CMAKE_MESSAGE_INDENT "  ")



# BLAS features
if(TARGET BLAS::BLAS)
  set(CMAKE_REQUIRED_LIBRARIES BLAS::BLAS)
else()
  set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
endif()
check_fortran_source_compiles("external :: zgemm3m; call zgemm3m(); end"
  BLAS_HAS_GEMM3M SRC_EXT F90)

unset(CMAKE_REQUIRED_LIBRARIES)



# LAPACK features
set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK)
check_fortran_source_compiles("external :: zheevr, dsyevr; call zheevr(); call dsyevr(); end"
  LAPACK_HAS_MRRR SRC_EXT F90)

# check for 2-stage solvers
set(_lapack_2stage_suffixes VX VD V)
if( LAPACK_HAS_MRRR )
  # add the MRRR 2-stage solver
  list(APPEND _lapack_2stage_suffixes VR)
endif()

set(_lapack_2stage TRUE)
foreach(suffix IN LISTS _lapack_2stage_suffixes)
  # check for names
  check_fortran_source_compiles("external :: zhee${suffix}_2stage
    external :: dsye${suffix}_2stage
    call zhee${suffix}_2stage()
    call dsye${suffix}_2stage()
    end" LAPACK_HAS_${suffix}_2STAGE SRC_EXT F90
    )
  if(NOT LAPACK_HAS_${suffix}_2STAGE)
    set(_lapack_2stage FALSE)
  endif()
endforeach()

if( ${_lapack_2stage} )
  set(LAPACK_HAS_2STAGE TRUE CACHE BOOL "2stage solvers available in LAPACK")
else()
  set(LAPACK_HAS_2STAGE FALSE CACHE BOOL "2stage solvers available in LAPACK")
endif()

unset(CMAKE_REQUIRED_LIBRARIES)


# Scalapack features
if(SIESTA_WITH_MPI)
  # Also check ScaLAPACK for MRRR

  # Odly enough adding LAPACK::LAPACK here will mess up the linking order
  # Since LAPACK::LAPACK is added in FindCustomScalapack as a dependency we do not needed.
  set(CMAKE_REQUIRED_LIBRARIES SCALAPACK::SCALAPACK)
  check_fortran_source_compiles("external :: pdsyevr, pzheevr; call pdsyevr(); call pzheevr(); end"
                               SCALAPACK_HAS_MRRR SRC_EXT F90)
  unset(CMAKE_REQUIRED_LIBRARIES)

  set(HAS_MRRR ${LAPACK_HAS_MRRR} AND ${SCALAPACK_HAS_MRRR} CACHE BOOL "MRRR routines available")

else()

  set(HAS_MRRR ${LAPACK_HAS_MRRR} CACHE BOOL "MRRR routines available")

endif(SIESTA_WITH_MPI)

list(POP_BACK CMAKE_MESSAGE_INDENT)
