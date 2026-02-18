set(top_srcdir "${PROJECT_SOURCE_DIR}/Src" )

set(sources

  atompla.f
  colinear.f
  denchar.f
  length.f
  local_reinit.f
  matvect.f
  neighb.f
  planed.f
  ranger.f
  readpla.f
  redata_denchar.f
  rhoofr.f
  timer_local.f
  wavofr.f
  wannofr.f
  wrout.f
)

list(
  APPEND
  sources

  ${top_srcdir}/precision.F
  ${top_srcdir}/bessph.f
  ${top_srcdir}/chkdim.f
  ${top_srcdir}/dismin.f
  ${top_srcdir}/dot.f
  ${top_srcdir}/m_io.f
  ${top_srcdir}/io.f
  ${top_srcdir}/alloc.F90
  ${top_srcdir}/files.f
  ${top_srcdir}/parallel.F
  ${top_srcdir}/printmatrix.F
  ${top_srcdir}/pxf.F90
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
  ${top_srcdir}/iodm.F
  ${top_srcdir}/memory_log.F90
  ${top_srcdir}/memory.F
  ${top_srcdir}/m_fft_gpfa.F
  ${top_srcdir}/radfft.f
  ${top_srcdir}/volcel.f
  ${top_srcdir}/debugmpi.F
  ${top_srcdir}/listsc.f
  ${top_srcdir}/atmparams.f
  ${top_srcdir}/atom_options.F90
  ${top_srcdir}/interpolation.f90
  ${top_srcdir}/xml.f
  ${top_srcdir}/radial.f
  ${top_srcdir}/atm_types.f
  ${top_srcdir}/spher_harm.f
  ${top_srcdir}/atmfuncs.f
  ${top_srcdir}/basis_types.f
  ${top_srcdir}/chemical.f
  ${top_srcdir}/basis_io.F
  ${top_srcdir}/qsort.F
  ${top_srcdir}/mmio.F
  ${top_srcdir}/pspltm1.F
  ${top_srcdir}/m_getopts.f90
  ${top_srcdir}/m_uuid.f90
  ${top_srcdir}/object_debug.F90
  ${top_srcdir}/class_Pair_Data1D.F90
  ${top_srcdir}/class_Fstack_Pair_Data1D.F90
)

siesta_add_executable(${PROJECT_NAME}.denchar
    ${sources}
  NAMESPACE_TARGET denchar
  )

target_link_libraries(
  ${PROJECT_NAME}.denchar
  PRIVATE
  libfdf::libfdf
  libpsml::libpsml
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
  )

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.denchar
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
