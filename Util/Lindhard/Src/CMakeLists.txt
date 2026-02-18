set(top_srcdir "${PROJECT_SOURCE_DIR}/Src" )

set(sources
  lindhard_kgrid.f90
  lindhard_reinit.f90
  lindhard_tools.f90
  lindhard.F90
)

list(
  APPEND
  sources
  ${top_srcdir}/alloc.F90
  ${top_srcdir}/cross.f
  ${top_srcdir}/cellsubs.f
  ${top_srcdir}/files.f
  ${top_srcdir}/idiag.f
  ${top_srcdir}/m_io.f
  ${top_srcdir}/minvec.f
  ${top_srcdir}/parallel.F
  ${top_srcdir}/precision.F
  ${top_srcdir}/reclat.f
  ${top_srcdir}/redcel.F
  ${top_srcdir}/sorting.f
  ${top_srcdir}/volcel.f
)


siesta_add_executable(${PROJECT_NAME}.lindhard
    ${sources}
  NAMESPACE_TARGET lindhard
  )
target_link_libraries(
  ${PROJECT_NAME}.lindhard
  PRIVATE
  $<$<BOOL:${SIESTA_WITH_OPENMP}>:OpenMP::OpenMP_Fortran>
  libfdf::libfdf
  xmlf90::xmlf90
  ${PROJECT_NAME}.libpsop
  ${PROJECT_NAME}.libncps
  libgridxc::libgridxc
  $<$<BOOL:${SIESTA_WITH_MPI}>:${PROJECT_NAME}.mpi>
  $<$<BOOL:${SIESTA_WITH_NETCDF}>:NetCDF::NetCDF_Fortran>
  $<$<BOOL:${SIESTA_WITH_NCDF}>:${PROJECT_NAME}.libncdf>
  ${PROJECT_NAME}.libfdict     # Even if not using NCDF
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
  $<$<BOOL:${SIESTA_NEEDS_ELPA_TARGET}>:Elpa::elpa>
  $<$<BOOL:${SIESTA_NEEDS_ELSI_TARGET_FOR_ELPA}>:elsi::elsi>
  $<$<BOOL:${SIESTA_WITH_MPI}>:SCALAPACK::SCALAPACK>
  LAPACK::LAPACK
)

# Ensure the version specification can find the auto-generated file
target_include_directories(
  ${PROJECT_NAME}.lindhard
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
  )

target_compile_definitions(
  ${PROJECT_NAME}.lindhard
  PRIVATE
  LINDHARD
  $<$<BOOL:${SIESTA_WITH_MPI}>:MPI>
  $<$<BOOL:${HAS_MRRR}>:SIESTA__MRRR>
  $<$<BOOL:${SIESTA_WITH_ELPA}>:SIESTA__ELPA>
  $<$<BOOL:${SIESTA_WITH_NETCDF}>:CDF>
  $<$<BOOL:${SIESTA_WITH_NCDF}>:NCDF NCDF_4>
  $<$<BOOL:${SIESTA_WITH_NCDF_PARALLEL}>:NCDF_PARALLEL>
)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.lindhard
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
