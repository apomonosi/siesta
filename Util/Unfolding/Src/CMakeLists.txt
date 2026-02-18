set(top_srcdir "${PROJECT_SOURCE_DIR}/Src" )

set(sources
  array.F90
  ${PROJECT_SOURCE_DIR}/Util/HSX/hsx_m.f90
  unfold.F90
  unfold_handlers_m.F90
)

list(
  APPEND
  sources

  ${top_srcdir}/precision.F
  ${top_srcdir}/broadcast_fdf_struct.F90
  ${top_srcdir}/m_io.f
  ${top_srcdir}/alloc.F90
  ${top_srcdir}/parallel.F
  ${top_srcdir}/memory_log.F90
  ${top_srcdir}/memory.F
  ${top_srcdir}/pxf.F90
  ${top_srcdir}/moreParallelSubs.F90
  ${top_srcdir}/m_mpi_utils.F
  ${top_srcdir}/interpolation.f90
  ${top_srcdir}/xml.f
  ${top_srcdir}/radial.f
  ${top_srcdir}/atm_types.f
  ${top_srcdir}/spher_harm.f
  ${top_srcdir}/atmfuncs.f
  ${top_srcdir}/reclat.f
  ${top_srcdir}/idiag.f
  ${top_srcdir}/cellsubs.f
  ${top_srcdir}/sorting.f
  ${top_srcdir}/minvec.f
  ${top_srcdir}/io.f
  ${top_srcdir}/files.f
  ${top_srcdir}/timestamp.f90
  ${top_srcdir}/m_walltime.f90
  ${top_srcdir}/m_wallclock.f90
  ${top_srcdir}/diag_option.F90
  ${PROJECT_BINARY_DIR}/Src/configured_values_m.F90
  ${top_srcdir}/diag.F90
  ${top_srcdir}/bloch_unfold.F90
  ${top_srcdir}/debugmpi.F
  ${top_srcdir}/find_kgrid.F
  ${top_srcdir}/posix_calls.f90
  ${top_srcdir}/atmparams.f
  ${top_srcdir}/basis_types.f
  ${top_srcdir}/atom_options.F90
  ${top_srcdir}/chemical.f
  ${top_srcdir}/basis_io.F
  ${top_srcdir}/siesta_geom.F90
  ${top_srcdir}/get_kpoints_scale.f90
  ${top_srcdir}/bessph.f
  ${top_srcdir}/m_fft_gpfa.F
  ${top_srcdir}/radfft.f
  ${top_srcdir}/periodic_table.f
  ${top_srcdir}/basis_specs.f
  ${top_srcdir}/hamann.f90
  ${top_srcdir}/m_filter.f90
  ${top_srcdir}/atom.F
  ${top_srcdir}/m_cite.F90
  ${top_srcdir}/dftu_specs.f
  ${top_srcdir}/broadcast_basis.F
  ${top_srcdir}/m_cell.f
  ${top_srcdir}/zmatrix.F
  ${top_srcdir}/coor.F
  ${top_srcdir}/chkdim.f
  ${top_srcdir}/volcel.f
  ${top_srcdir}/digcel.f
  ${top_srcdir}/arw.f
  ${top_srcdir}/redcel.F
  ${top_srcdir}/printmatrix.F
  ${top_srcdir}/schecomm.F
  ${top_srcdir}/class_Geometry.F90
  ${top_srcdir}/class_Data2D.F90
  ${top_srcdir}/class_OrbitalDistribution.F90
  ${top_srcdir}/class_Sparsity.F90
  ${top_srcdir}/class_SpData2D.F90
  ${top_srcdir}/class_Pair_Geometry_SpData2D.F90
  ${top_srcdir}/class_Fstack_Pair_Geometry_SpData2D.F90
  ${top_srcdir}/class_Data1D.F90
  ${top_srcdir}/class_SpData1D.F90
  ${top_srcdir}/sparse_matrices.F90
  ${top_srcdir}/domain_decom.F
  ${top_srcdir}/spatial.F
  ${top_srcdir}/parallelsubs.F
  ${top_srcdir}/pspltm1.F
  ${top_srcdir}/mmio.F
  ${top_srcdir}/qsort.F
  ${top_srcdir}/timer_tree.f90
  ${top_srcdir}/m_timer.F90
  ${top_srcdir}/timer.F90
  ${top_srcdir}/m_uuid.f90
  ${top_srcdir}/object_debug.F90
  ${top_srcdir}/nag.f
  ${top_srcdir}/m_spin.F90
  ${top_srcdir}/m_vee_integrals.F90
)

siesta_add_executable(${PROJECT_NAME}.unfold
    ${sources}
  NAMESPACE_TARGET unfold
  )

# This seems to be necessary for Ninja, but only
# when SIESTA_NEEDS_ELSI_TARGET_FOR_ELPA is true...

target_include_directories(
  ${PROJECT_NAME}.unfold
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
  )

# Putting the MPI dependency first seems to fix some issues (see Src/CMakeLists).
target_link_libraries(
  ${PROJECT_NAME}.unfold
  PRIVATE
  $<$<BOOL:${SIESTA_WITH_MPI}>:${PROJECT_NAME}.mpi>
  $<$<BOOL:${SIESTA_WITH_NETCDF}>:NetCDF::NetCDF_Fortran>
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libpsop
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
  libfdf::libfdf
  libpsml::libpsml
  libgridxc::libgridxc
  $<$<BOOL:${SIESTA_NEEDS_ELPA_TARGET}>:Elpa::elpa>
  $<$<BOOL:${SIESTA_NEEDS_ELSI_TARGET_FOR_ELPA}>:elsi::elsi>
  $<$<BOOL:${SIESTA_WITH_MPI}>:SCALAPACK::SCALAPACK>
  LAPACK::LAPACK
  $<$<BOOL:${SIESTA_WITH_OPENMP}>:OpenMP::OpenMP_Fortran>
 )

target_compile_definitions(
  ${PROJECT_NAME}.unfold
  PRIVATE
  $<$<BOOL:${SIESTA_WITH_MPI}>:MPI>
  $<$<BOOL:${HAS_MRRR}>:SIESTA__MRRR>
  $<$<BOOL:${LAPACK_HAS_2STAGE}>:SIESTA__DIAG_2STAGE>
  $<$<BOOL:${SIESTA_WITH_ELPA}>:SIESTA__ELPA>
  $<$<BOOL:${SIESTA_WITH_NETCDF}>:CDF>
)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.unfold
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
endif()
