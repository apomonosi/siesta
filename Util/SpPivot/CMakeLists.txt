#
set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

siesta_add_executable(${PROJECT_NAME}.pvtsp
  NAMESPACE_TARGET pvtsp

   # Required (but not used) dependencies
   ${top_srcdir}/atm_types.f
   ${top_srcdir}/atmfuncs.f
   ${top_srcdir}/atmparams.f
   ${top_srcdir}/atom.F
   ${top_srcdir}/atomlist.f
   ${top_srcdir}/atom_options.F90
   ${top_srcdir}/basis_specs.f
   ${top_srcdir}/basis_types.f
   ${top_srcdir}/bessph.f
   ${top_srcdir}/cellsubs.f
   ${top_srcdir}/chkdim.f
   ${top_srcdir}/class_Geometry.F90
   ${top_srcdir}/class_Pair_Geometry_SpData2D.F90
   ${top_srcdir}/class_Fstack_Pair_Geometry_SpData2D.F90
   ${top_srcdir}/chemical.f
   ${top_srcdir}/dftu_specs.f
   ${top_srcdir}/domain_decom.F
   ${top_srcdir}/files.f
   ${top_srcdir}/find_kgrid.F
   ${top_srcdir}/hamann.f90
   ${top_srcdir}/idiag.f
   ${top_srcdir}/interpolation.f90
   ${top_srcdir}/kpoint_convert.f90
   ${top_srcdir}/kpoint_t.F90
   ${top_srcdir}/minvec.f
   ${top_srcdir}/mmio.F
   ${top_srcdir}/m_cite.F90
   ${top_srcdir}/m_char.f90
   ${top_srcdir}/m_fft_gpfa.F
   ${top_srcdir}/m_filter.f90
   ${top_srcdir}/m_spin.F90
   ${top_srcdir}/m_vee_integrals.F90
   ${top_srcdir}/parallelsubs.F
   ${top_srcdir}/periodic_table.f
   ${top_srcdir}/printmatrix.F
   ${top_srcdir}/pspltm1.F
   ${top_srcdir}/qsort.F
   ${top_srcdir}/siesta_geom.F90
   ${top_srcdir}/sorting.f
   ${top_srcdir}/sparse_matrices.F90
   ${top_srcdir}/spatial.F
   ${top_srcdir}/radial.f
   ${top_srcdir}/radfft.f
   ${top_srcdir}/spher_harm.f
   ${top_srcdir}/siesta_cml.f90
   ${top_srcdir}/schecomm.F
   ${top_srcdir}/xml.f

   ${top_srcdir}/precision.F
   ${top_srcdir}/parallel.F
   ${top_srcdir}/alloc.F90
   ${top_srcdir}/m_io.f
   ${top_srcdir}/reclat.f
   ${top_srcdir}/io.f
   ${top_srcdir}/memory_log.F90
   ${top_srcdir}/memory.F
   ${top_srcdir}/m_uuid.f90
   ${top_srcdir}/object_debug.F90
   ${top_srcdir}/class_Sparsity.F90
   ${top_srcdir}/class_OrbitalDistribution.F90
   ${top_srcdir}/class_Data1D.F90
   ${top_srcdir}/class_Data2D.F90
   ${top_srcdir}/class_SpData1D.F90
   ${top_srcdir}/class_SpData2D.F90
   ${top_srcdir}/intrinsic_missing.F90
   ${top_srcdir}/geom_helper.f90
   ${top_srcdir}/io_sparse.F90
   ${top_srcdir}/m_os.F90
   ${top_srcdir}/create_Sparsity_SC.F90
   ${top_srcdir}/m_sparse.F90
   ${top_srcdir}/io_hsx.F90
   ${top_srcdir}/ts_io_hs.F90
   ${top_srcdir}/m_region.F90
   ${top_srcdir}/m_pivot_methods.F90
   ${top_srcdir}/m_pivot.F90
   ${top_srcdir}/m_interpolate.F90
   ${top_srcdir}/m_sparsity_handling.F90

   pvtsp.F90
   local_timer.f90
)

target_link_libraries(${PROJECT_NAME}.pvtsp PRIVATE ${PROJECT_NAME}.libsys
    libfdf::libfdf
    libpsml::libpsml
    libgridxc::libgridxc
    ${PROJECT_NAME}.libncps
    ${PROJECT_NAME}.libpsop
    ${PROJECT_NAME}.libunits
    ${PROJECT_NAME}.mpi_macros
    LAPACK::LAPACK
)


if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.pvtsp
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

