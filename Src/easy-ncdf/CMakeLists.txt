#
# Use a pre-generated source file (see top of file)
#
siesta_add_library(
  ${PROJECT_NAME}.libncdf
  NO_LINKER_FLAGS
  STATIC 
    netcdf_ncdf.F90
  )

if(SIESTA_WITH_NCDF_PARALLEL)
  target_compile_definitions(
    ${PROJECT_NAME}.libncdf
    PRIVATE
    NCDF_PARALLEL
  )

  # Note that libncdf uses the MPI module directly, through 'use mpi'
  # This dependency has been put first to avoid problems with some
  # MPI installations (e.g. when using oneapi-mpi with the GNU compiler in CINECA)

  target_link_libraries(
     ${PROJECT_NAME}.libncdf
     PRIVATE
     MPI::MPI_Fortran
  )

endif()

target_link_libraries(
   ${PROJECT_NAME}.libncdf
   # 'PUBLIC' to bring in fdict
   PRIVATE
   ${PROJECT_NAME}.libfdict
   NetCDF::NetCDF_Fortran
)

target_include_directories(
  ${PROJECT_NAME}.libncdf
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

