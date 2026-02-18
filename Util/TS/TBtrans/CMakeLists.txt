set(top_srcdir "${PROJECT_SOURCE_DIR}/Src" )

set(sources

 tbt_handlers_m.F90
 m_tbt_contour.F90
 m_tbt_dH.F90
 m_tbt_dSE.F90
 m_tbt_delta.F90
 m_tbt_diag.F90
 m_tbt_gf.F90
 m_tbt_hs.F90
 m_tbt_kpoint.F90
 m_tbt_kregions.F90
 m_tbt_options.F90
 m_tbt_proj.F90
 m_tbt_regions.F90
 m_tbt_save.F90
 m_tbt_sigma_save.F90
 m_tbt_sparse_helper.F90
 m_tbt_tri_init.F90
 m_tbt_tri_scat.F90
 m_tbt_trik.F90
 m_tbtrans.F90
 m_verbosity.f90
 sorted_search.F90
 tbt_end.F90
 tbt_init.F90
 tbt_init_output.f90
 tbt_reinit_m.F90
 tbtrans.F90

)

list(
  APPEND
  sources

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
  ${top_srcdir}/chkdim.f
  ${top_srcdir}/class_Geometry.F90
  ${top_srcdir}/class_Pair_Geometry_SpData2D.F90
  ${top_srcdir}/class_Fstack_Pair_Geometry_SpData2D.F90
  ${top_srcdir}/chemical.f
  ${top_srcdir}/dftu_specs.f
  ${top_srcdir}/domain_decom.F
  ${top_srcdir}/hamann.f90
  ${top_srcdir}/interpolation.f90
  ${top_srcdir}/mmio.F
  ${top_srcdir}/m_fft_gpfa.F
  ${top_srcdir}/m_filter.f90
  ${top_srcdir}/m_vee_integrals.F90
  ${top_srcdir}/parallelsubs.F
  ${top_srcdir}/periodic_table.f
  ${top_srcdir}/printmatrix.F
  ${top_srcdir}/pspltm1.F
  ${top_srcdir}/qsort.F
  ${top_srcdir}/sparse_matrices.F90
  ${top_srcdir}/spatial.F
  ${top_srcdir}/radial.f
  ${top_srcdir}/radfft.f
  ${top_srcdir}/spher_harm.f
  ${top_srcdir}/schecomm.F
  ${top_srcdir}/xml.f

  ${top_srcdir}/alloc.F90
  ${top_srcdir}/basic_func.inc
  ${top_srcdir}/basic_type.inc
  ${top_srcdir}/bloch_unfold.F90
  ${top_srcdir}/broadcast_fdf_struct.F90
  ${top_srcdir}/byte_count.F90
  ${top_srcdir}/cellsubs.f
  ${top_srcdir}/class_Data1D.F90
  ${top_srcdir}/class_Data1D.T90
  ${top_srcdir}/class_Data2D.F90
  ${top_srcdir}/class_Data2D.T90
  ${top_srcdir}/class_Distribution.F90
  ${top_srcdir}/class_OrbitalDistribution.F90
  ${top_srcdir}/class_SpData1D.F90
  ${top_srcdir}/class_SpData1D.T90
  ${top_srcdir}/class_SpData2D.F90
  ${top_srcdir}/class_SpData2D.T90
  ${top_srcdir}/class_Sparsity.F90
  ${top_srcdir}/class_TriMat.F90
  ${top_srcdir}/class_TriMat.T90
  ${top_srcdir}/cli_m.f90
  ${PROJECT_BINARY_DIR}/Src/configured_values_m.F90
  ${top_srcdir}/create_Sparsity_SC.F90
  ${top_srcdir}/create_Sparsity_Union.F90
  ${top_srcdir}/cross.f
  ${top_srcdir}/debugmpi.F
  ${top_srcdir}/densematrix.f90
  ${top_srcdir}/diag.F90
  ${top_srcdir}/diag_option.F90
  ${top_srcdir}/fdf_extra.F90
  ${top_srcdir}/files.f
  ${top_srcdir}/find_kgrid.F
  ${top_srcdir}/geom_helper.f90
  ${top_srcdir}/get_kpoints_scale.f90
  ${top_srcdir}/idiag.f
  ${top_srcdir}/intrinsic_missing.F90
  ${top_srcdir}/io.f
  ${top_srcdir}/kpoint_convert.f90
  ${top_srcdir}/kpoint_t.F90
  ${top_srcdir}/m_char.f90
  ${top_srcdir}/m_cite.F90
  ${top_srcdir}/m_gauss_fermi_17.f90
  ${top_srcdir}/m_gauss_fermi_18.f90
  ${top_srcdir}/m_gauss_fermi_19.f90
  ${top_srcdir}/m_gauss_fermi_20.f90
  ${top_srcdir}/m_gauss_fermi_22.f90
  ${top_srcdir}/m_gauss_fermi_24.f90
  ${top_srcdir}/m_gauss_fermi_26.f90
  ${top_srcdir}/m_gauss_fermi_28.f90
  ${top_srcdir}/m_gauss_fermi_30.f90
  ${top_srcdir}/m_gauss_fermi_inf.f90
  ${top_srcdir}/m_gauss_quad.f90
  ${top_srcdir}/m_geom_aux.f90
  ${top_srcdir}/m_geom_box.f90
  ${top_srcdir}/m_geom_coord.f90
  ${top_srcdir}/m_geom_objects.f90
  ${top_srcdir}/m_geom_plane.f90
  ${top_srcdir}/m_geom_square.f90
  ${top_srcdir}/m_handle_sparse.F90
  ${top_srcdir}/m_integrate.f90
  ${top_srcdir}/m_interpolate.F90
  ${top_srcdir}/m_io.f
  ${top_srcdir}/io_sparse.F90
  ${top_srcdir}/m_iodm.F90
  ${top_srcdir}/m_iterator.f90
  ${top_srcdir}/m_mat_invert.F90
  ${top_srcdir}/m_mpi_utils.F
  ${top_srcdir}/m_os.F90
  ${top_srcdir}/m_pivot.F90
  ${top_srcdir}/m_pivot_array.f90
  ${top_srcdir}/m_pivot_methods.F90
  ${top_srcdir}/m_region.F90
  ${top_srcdir}/m_sparse.F90
  ${top_srcdir}/m_sparsity_handling.F90
  ${top_srcdir}/m_spin.F90
  ${top_srcdir}/m_timer.F90
  ${top_srcdir}/m_trimat_invert.F90
  ${top_srcdir}/m_ts_aux.F90
  ${top_srcdir}/m_ts_cctype.f90
  ${top_srcdir}/m_ts_chem_pot.F90
  ${top_srcdir}/m_ts_contour.f90
  ${top_srcdir}/m_ts_contour_eq.f90
  ${top_srcdir}/m_ts_contour_neq.f90
  ${top_srcdir}/m_ts_debug.F90
  ${top_srcdir}/m_ts_elec_se.F90
  ${top_srcdir}/m_ts_electrode.F90
  ${top_srcdir}/m_ts_method.f90
  ${top_srcdir}/m_ts_gf.F90
  ${top_srcdir}/ts_io_hs.F90
  ${top_srcdir}/io_hsx.F90
  ${top_srcdir}/m_ts_io_contour.f90
  ${top_srcdir}/m_ts_io_ctype.f90
  ${top_srcdir}/m_ts_iodm.F90
  ${top_srcdir}/m_ts_pivot.F90
  ${top_srcdir}/m_ts_rgn2trimat.F90
  ${top_srcdir}/m_ts_sparse.F90
  ${top_srcdir}/m_ts_sparse_helper.F90
  ${top_srcdir}/m_ts_tdir.f90
  ${top_srcdir}/m_ts_tri_common.F90
  ${top_srcdir}/m_ts_tri_scat.F90
  ${top_srcdir}/m_ts_trimat_invert.F90
  ${top_srcdir}/m_uuid.f90
  ${top_srcdir}/m_wallclock.f90
  ${top_srcdir}/m_walltime.f90
  ${top_srcdir}/memory.F
  ${top_srcdir}/memory_log.F90
  ${top_srcdir}/minvec.f
  ${top_srcdir}/moreParallelSubs.F90
  ${top_srcdir}/ncdf_io.F90
  ${top_srcdir}/object_debug.F90
  ${top_srcdir}/parallel.F
  ${top_srcdir}/posix_calls.f90
  ${top_srcdir}/precision.F
  ${top_srcdir}/pxf.F90
  ${top_srcdir}/reclat.f
  ${top_srcdir}/runinfo_m.F90
  ${top_srcdir}/siesta_geom.F90
  ${top_srcdir}/siesta_version_info.F90
  ${top_srcdir}/siesta_cml.f90
  ${top_srcdir}/sorting.f
  ${top_srcdir}/spin_subs.F90
  ${top_srcdir}/timer.F90
  ${top_srcdir}/timer_tree.f90
  ${top_srcdir}/timestamp.f90
  ${top_srcdir}/ts_electrode.F90
)


siesta_add_executable(${PROJECT_NAME}.tbtrans
    ${sources}
  NAMESPACE_TARGET tbtrans
  )
target_link_libraries(
  ${PROJECT_NAME}.tbtrans
  PRIVATE
  $<$<BOOL:${SIESTA_WITH_OPENMP}>:OpenMP::OpenMP_Fortran>
  libfdf::libfdf
  xmlf90::xmlf90 # required dependency, not used by TBTRANS!
  ${PROJECT_NAME}.libpsop # required dependency, not used by TBTRANS!
  ${PROJECT_NAME}.libncps # required dependency, not used by TBTRANS!
  libgridxc::libgridxc # required dependency, not used by TBTRANS!
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
  ${PROJECT_NAME}.tbtrans
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
  )

target_compile_definitions(
  ${PROJECT_NAME}.tbtrans
  PRIVATE
  TBTRANS # Notify about TBTRANS
  $<$<BOOL:${SIESTA_WITH_MPI}>:MPI>
  $<$<BOOL:${HAS_MRRR}>:SIESTA__MRRR>
  $<$<BOOL:${SIESTA_WITH_ELPA}>:SIESTA__ELPA>
  $<$<BOOL:${SIESTA_WITH_NETCDF}>:CDF>
  $<$<BOOL:${SIESTA_WITH_NCDF}>:NCDF NCDF_4>
  $<$<BOOL:${SIESTA_WITH_NCDF_PARALLEL}>:NCDF_PARALLEL>
)

if( CMAKE_Fortran_COMPILER_ID MATCHES NVHPC )

  set_property(
    SOURCE
    ${top_srcdir}/m_ts_electrode.F90
    APPEND PROPERTY COMPILE_OPTIONS
    "-Mnoinline"
  )
endif()


# Create string with line breaks and continuations in order to not exceed
# the maximum line length (132 characters) of Fortran 2003.
# There is not check that the maximum number of continuation lines
# (255 in Fortran 2003) is not exceeded.
siesta_conda_fix(Fortran_FLAGS_CURRENT)
siesta_get_multiline(LENGTH 128
  Fortran_FLAGS_CURRENT
  OUTPUT Fortran_FLAGS_CURRENT_multiline
  )

# This isn't used since it may contain generator
# expressions, and the multiline method can cut it
# at bad positions.
get_target_property(
  Fortran_PPFLAGS_CURRENT
  ${PROJECT_NAME}.tbtrans
  COMPILE_DEFINITIONS
  )

siesta_conda_fix(Fortran_PPFLAGS_CURRENT)
siesta_get_multiline(LENGTH 128
  Fortran_PPFLAGS_CURRENT
  OUTPUT Fortran_PPFLAGS_CURRENT_multiline
  )

# Create version file with the version string implemented
# TODO this will only be generated at cmake generation, i.e.
#  git pull | cmake --build_build -- siesta-siesta
# will not regenerate the version information.
configure_file(
  "${PROJECT_SOURCE_DIR}/Src/version-info-template.inc"
  version-info.inc
  @ONLY
  )


if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.tbtrans
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
