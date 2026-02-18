#
# This function uses 'SIESTA_TESTS_MPI_NUMPROCS' as set in a previous
# call to TestMPIConfig.
#

# Note TARGET_FILE idiom for generality. This also covers the case in
# which the current directory is not in the PATH ('./prog' could work
# too). Otherwise mpiexec would not find it.

function(setup_test_mpi prog)

   set(mpiexec_args
          ${MPIEXEC_NUMPROC_FLAG} ${SIESTA_TESTS_MPI_NUMPROCS}
          ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:${prog}>
          ${MPIEXEC_POSTFLAGS})

   add_test( NAME    ${prog}_mpi_np_${SIESTA_TESTS_MPI_NUMPROCS}
             COMMAND ${MPIEXEC_EXECUTABLE} ${mpiexec_args} )

endfunction()

siesta_add_executable(${PROJECT_NAME}.pi3
  pi3.F
  )

target_link_libraries(${PROJECT_NAME}.pi3
  PRIVATE
    ${PROJECT_NAME}.mpi
    ${PROJECT_NAME}.mpi_macros
)
#
# To check "singleton" MPI runs (i.e., that we can run
# a parallel program without 'mpirun'). This is not
# possible in some systems (e.g. IBM's MPI)
#
add_test(NAME pi3_Runs_Singleton COMMAND ${PROJECT_NAME}.pi3)

# Normal MPI run with ${SIESTA_TESTS_MPI_NUMPROCS}
setup_test_mpi(${PROJECT_NAME}.pi3)


siesta_add_executable(${PROJECT_NAME}.newcomm
  newcomm.F
  ${PROJECT_SOURCE_DIR}/Src/MPI
)

target_link_libraries(${PROJECT_NAME}.newcomm
  PRIVATE
    ${PROJECT_NAME}.mpi
    ${PROJECT_NAME}.mpi_macros
  )
setup_test_mpi(${PROJECT_NAME}.newcomm)

#
siesta_add_executable(${PROJECT_NAME}.blacs_prb blacs_prb.f90 )
target_link_libraries(${PROJECT_NAME}.blacs_prb
  PRIVATE
  SCALAPACK::SCALAPACK
  MPI::MPI_Fortran
  LAPACK::LAPACK
)
setup_test_mpi(${PROJECT_NAME}.blacs_prb)


siesta_add_executable(${PROJECT_NAME}.pblas_prb  pblas_prb.f90 )
target_link_libraries(${PROJECT_NAME}.pblas_prb
  PRIVATE
  SCALAPACK::SCALAPACK
  MPI::MPI_Fortran
  LAPACK::LAPACK
)

file(COPY "${CMAKE_CURRENT_LIST_DIR}/pblas.dat" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

#
# Run the test automatically only if we have enough ranks available.
#
if("${SIESTA_TESTS_MPI_NUMPROCS}" GREATER_EQUAL 4)
  setup_test_mpi(${PROJECT_NAME}.pblas_prb)
endif()


