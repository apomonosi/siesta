set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

set(sources
  assignation.F90
  centerall.F90
  elecfield.F90
  coulomb.F90
  functions.F90
  graphite.f90
  harmonic_restraints.f90
  ioxv_restr.F90
  libsiesta_init.F90
  linkatoms.F90
  ljef.F90
  mm_assign.f90
  mm_ene_fce.f90
  mm_topology.f90
  mm_units.f90
  pbc.f90
  qmmm_files.f90
  qmmm_fixed.F90
  qmmm_fsiesta.F90
  qmmm_ioxv.f90
  qmmm_lst_blk.f90
  qmmm_mneighb.f90
  qmmm_reinit.F90
  qmmm_siesta_move.F90
  qmmm_struct_init.F90
  qmmm_timer.f90
  readall.F90
  reinsert_atoms.f90
  siesta_qmmm.F90
  siesta_qmmm_options.F90
  version.F90
  writeall.F90

  replaced_siesta_modules/atomlist.f90
  replaced_siesta_modules/atm_types.f90
  replaced_siesta_modules/atmfuncs.f90
  replaced_siesta_modules/empty_modules.f90
)

list(
  APPEND
  sources
  # Actually, few of the modules here are actually required, the majority
  # hop-on due to dependencies within those same files.

  # The following files are required for the QMMM driver to use the same
  # molecular mechanics routines as the rest of SIESTA.
  ${top_srcdir}/broyden_optim.F
  ${top_srcdir}/cell_broyden_optim.F
  ${top_srcdir}/cell_fire_optim.F
  ${top_srcdir}/cellsubs.f
  ${top_srcdir}/cgvc.F
  ${top_srcdir}/cgvc_zmatrix.F
  ${top_srcdir}/conjgr.f
  ${top_srcdir}/conjgr_old.f
  ${top_srcdir}/dynamics.F
  ${top_srcdir}/fire_optim.F
  ${top_srcdir}/m_fire.f90
  ${top_srcdir}/geom_helper.f90
  ${top_srcdir}/iomd.f
  ${top_srcdir}/iotdxv.F
  ${top_srcdir}/ioxv.F
  ${top_srcdir}/iozm.F
  ${top_srcdir}/m_fixed.F90
  ${top_srcdir}/m_broyddj.f90
  ${top_srcdir}/m_broyddj_nocomm.f90
  ${top_srcdir}/metaforce.F
  ${top_srcdir}/ofc.f90
  ${top_srcdir}/zm_broyden_optim.F
  ${top_srcdir}/zm_fire_optim.F
  ${top_srcdir}/zmatrix.F
  ${top_srcdir}/vmb.F
  ${top_srcdir}/constr.f

  # These are used for communications.
  ${top_srcdir}/fsockets.f90


  # File IO
  ${top_srcdir}/broadcast_fdf_struct.F90
  ${top_srcdir}/fdf_extra.F90
  ${top_srcdir}/files.f
  ${top_srcdir}/io.f
  ${top_srcdir}/iocg.f
  ${top_srcdir}/iofa.f90
  ${top_srcdir}/m_io.f
  ${top_srcdir}/m_iostruct.f
  ${top_srcdir}/md_out.F90
  ${top_srcdir}/m_uuid.f90

  # Mainly used for variable storage.
  ${top_srcdir}/chemical.f
  ${top_srcdir}/m_energies.F90
  ${top_srcdir}/m_forces.F90
  ${top_srcdir}/m_kinetic.F90
  ${top_srcdir}/m_steps.F90
  ${top_srcdir}/m_stress.F90
  ${top_srcdir}/m_target_stress.F90

  # Utilitarian functions
  ${top_srcdir}/alloc.F90
  ${top_srcdir}/alloc_handlers_m.F90
  ${top_srcdir}/chkdim.f
  ${top_srcdir}/cli_m.f90
  ${top_srcdir}/m_check_walltime.f90
  ${top_srcdir}/m_cite.F90
  ${top_srcdir}/m_timer.F90
  ${top_srcdir}/memory.F
  ${top_srcdir}/memory_all.F90
  ${top_srcdir}/memory_log.F90
  ${top_srcdir}/memory_snapshot.f90
  ${top_srcdir}/m_wallclock.f90
  ${top_srcdir}/m_walltime.f90


  # These are miscellaneous dependencies due to the use of others.
  # Mainly, read_options.F90. For example, read_options requires m_charge_add.F90,
  # which in turns requires all of m_geom_*.f90 files. Or, in another case,
  # m_fixed.f90 requires m_region.F90, which itself requires class_sparsity.F90
  # and class_OrbitalDistribution.F90.
  ${top_srcdir}/chemical.f
  ${top_srcdir}/class_OrbitalDistribution.F90
  ${top_srcdir}/class_Sparsity.F90

  ${top_srcdir}/densematrix.f90
  ${top_srcdir}/digcel.f
  ${top_srcdir}/idiag.f
  ${top_srcdir}/intrinsic_missing.F90
  ${top_srcdir}/inver.f

  ${top_srcdir}/m_cell.f
  ${top_srcdir}/m_char.f90
  ${top_srcdir}/m_charge_add.F90

  ${top_srcdir}/m_geom_aux.f90
  ${top_srcdir}/m_geom_box.f90
  ${top_srcdir}/m_geom_coord.f90
  ${top_srcdir}/m_geom_objects.f90
  ${top_srcdir}/m_geom_plane.f90
  ${top_srcdir}/m_geom_square.f90

  ${top_srcdir}/m_mesh_node.F90
  ${top_srcdir}/m_mpi_utils.F
  ${top_srcdir}/m_region.F90
  ${top_srcdir}/m_spin.F90

  ${top_srcdir}/m_ts_global_vars.f90
  ${top_srcdir}/madelung.f

  ${top_srcdir}/molecularmechanics.F90
  ${top_srcdir}/moreParallelSubs.F90
  ${top_srcdir}/parallel.F
  ${top_srcdir}/periodic_table.f
  ${top_srcdir}/pixmol.f
  ${top_srcdir}/precision.F
  ${top_srcdir}/pxf.F90
  ${top_srcdir}/randomg.f90
  ${top_srcdir}/reclat.f
  ${top_srcdir}/redcel.F
  ${top_srcdir}/reinit_m.F90
  ${top_srcdir}/remove_intramol_pressure.f90
  ${top_srcdir}/reord.f
  ${top_srcdir}/rusage.f90
  ${top_srcdir}/siesta_cml.f90
  ${top_srcdir}/siesta_cmlsubs.F90
  ${top_srcdir}/siesta_geom.F90
  ${top_srcdir}/siesta_handlers_m.F90
  ${top_srcdir}/siesta_options.F90
  ${top_srcdir}/siesta_version_info.F90
  ${top_srcdir}/sockets.c
  ${top_srcdir}/timer.F90
  ${top_srcdir}/timer_tree.f90
  ${top_srcdir}/timestamp.f90
  ${top_srcdir}/typecell.f
  ${top_srcdir}/uncell.f
  ${top_srcdir}/volcel.f
  ${top_srcdir}/write_md_record.F
  ${top_srcdir}/xml.f
)

siesta_add_executable(
  ${PROJECT_NAME}.siesta_qmmm
  ${sources}
)

install(
  TARGETS ${PROJECT_NAME}.siesta_qmmm
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

configure_file(
  ${top_srcdir}/version-info-template.inc
  ${CMAKE_CURRENT_BINARY_DIR}/generated/version-info.inc
  COPYONLY
)

target_include_directories(
  ${PROJECT_NAME}.siesta_qmmm
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/generated
)

target_compile_definitions(
    ${PROJECT_NAME}.siesta_qmmm
    PUBLIC
    "$<$<BOOL:${SIESTA_WITH_MPI}>:MPI>"
    "$<$<BOOL:${SIESTA_WITH_GRID_SP}>:GRID_SP>"
    "$<$<BOOL:${HAS_MRRR}>:SIESTA__MRRR>"
    "$<$<BOOL:${SIESTA_WITH_ELPA}>:SIESTA__ELPA>"
    "$<$<BOOL:${SIESTA_WITH_NETCDF}>:CDF>"
    "$<$<BOOL:${SIESTA_WITH_NCDF}>:NCDF>"
    "$<$<BOOL:${SIESTA_WITH_NCDF}>:NCDF_4>"
    "$<$<BOOL:${SIESTA_WITH_NCDF_PARALLEL}>:NCDF_PARALLEL>"
    "$<$<BOOL:${SIESTA_WITH_LUA}>:SIESTA__FLOOK>"
    "$<$<BOOL:${SIESTA_WITH_DFTD3}>:SIESTA__DFTD3>"
)

target_link_libraries(
  ${PROJECT_NAME}.siesta_qmmm
  PRIVATE
  $<$<BOOL:${SIESTA_WITH_MPI}>:${PROJECT_NAME}.mpi>
  libfdf::libfdf
  libgridxc::libgridxc
  libpsml::libpsml
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libpsop
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
  $<$<BOOL:${SIESTA_WITH_FLOOK}>:flook::flook>
  $<$<OR:$<BOOL:${SIESTA_WITH_NCDF}>,$<BOOL:${SIESTA_WITH_FLOOK}>>:${PROJECT_NAME}.libfdict>
  $<$<BOOL:${SIESTA_WITH_NCDF}>:${PROJECT_NAME}.libncdf>
  $<$<BOOL:${SIESTA_WITH_NETCDF}>:NetCDF::NetCDF_Fortran>
  $<$<BOOL:${SIESTA_WITH_DFTD3}>:s-dftd3::s-dftd3>
  LAPACK::LAPACK
  $<$<BOOL:${SIESTA_WITH_OPENMP}>:OpenMP::OpenMP_Fortran>
 )


