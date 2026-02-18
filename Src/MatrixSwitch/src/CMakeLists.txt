# Built-in MatrixSwitch
#
siesta_add_library(${PROJECT_NAME}.libms
  STATIC

  MatrixSwitch.F90
  MatrixSwitch_m_add.F90
  MatrixSwitch_m_copy.F90
  MatrixSwitch_m_register.F90
  MatrixSwitch_m_set.F90
  MatrixSwitch_mm_multiply.F90
  MatrixSwitch_ops.F90
  MatrixSwitch_wrapper.F90
  MatrixSwitch_wrapper_params.F90
)

if (SIESTA_WITH_MPI)

  # specific preprocessor variables for compilation of MS itself
  target_compile_definitions(
   ${PROJECT_NAME}.libms
   PRIVATE
   HAVE_MPI
   HAVE_SCALAPACK
  )

  # Note that this library does NOT use mpi_siesta... We need the
  # standard MPI dependency, for itself and any users

  target_link_libraries(${PROJECT_NAME}.libms
    PUBLIC
      SCALAPACK::SCALAPACK
  )

# Note that the variables below are already incorporated as
# properties of the MPI::MPI_Fortran target by the FindMPI module. For
# example:
# set_property(TARGET MPI::MPI_${LANG}
#                     PROPERTY INTERFACE_COMPILE_DEFINITIONS
#                     "${MPI_${LANG}_COMPILE_DEFINITIONS}")

## --- old ---
## This enables correct compilation and also correct flags for
## MPI builds.
## Also, it ensures that the include directories are correctly setup
##  target_compile_options(${PROJECT_NAME}.libms PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
##  target_include_directories(${PROJECT_NAME}.libms PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
##  target_link_directories(${PROJECT_NAME}.libms PUBLIC ${MPI_Fortran_LIBRARY_DIRS})
##  target_link_libraries(${PROJECT_NAME}.libms PUBLIC  ${MPI_Fortran_LIBRARIES})

else()

  # specific preprocessor variables for compilation of MS itself
  target_compile_definitions(
   ${PROJECT_NAME}.libms
   PRIVATE
   HAVE_LAPACK
  )

  # add LAPACK dependency
  target_link_libraries(${PROJECT_NAME}.libms
    PUBLIC
      LAPACK::LAPACK
  )

endif(SIESTA_WITH_MPI)

# So that targets using MS can find its modules:
target_include_directories(${PROJECT_NAME}.libms
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}"
)

