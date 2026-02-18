set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

set(sources)

list(
    APPEND
    sources

    "${top_srcdir}/alloc.F90"
    "${top_srcdir}/arw.f"
    "${top_srcdir}/atm_types.f"
    "${top_srcdir}/atmparams.f"
    "${top_srcdir}/atom.F"
    "${top_srcdir}/atom_options.F90"
    "${top_srcdir}/basis_io.F"
    "${top_srcdir}/basis_specs.f"
    "${top_srcdir}/basis_types.f"
    "${top_srcdir}/bessph.f"
    "${top_srcdir}/cellsubs.f"
    "${top_srcdir}/chemical.f"
    "${top_srcdir}/dftu_specs.f"
    "${top_srcdir}/dot.f"
    "${top_srcdir}/files.f"
    "${top_srcdir}/get_kpoints_scale.f90"
    "${top_srcdir}/hamann.f90"
    "${top_srcdir}/interpolation.f90"
    "${top_srcdir}/io.f"
    "${top_srcdir}/minvec.f"
    "${top_srcdir}/m_cite.F90"
    "${top_srcdir}/m_fft_gpfa.F"
    "${top_srcdir}/m_filter.f90"
    "${top_srcdir}/m_io.f"
    "${top_srcdir}/m_spin.F90"
    "${top_srcdir}/m_walltime.f90"
    "${top_srcdir}/memory.F"
    "${top_srcdir}/memory_log.F90"
    "${top_srcdir}/nag.f"
    "${top_srcdir}/parallel.F"
    "${top_srcdir}/parsing.f"
    "${top_srcdir}/periodic_table.f"
    "${top_srcdir}/precision.F"
    "${top_srcdir}/pxf.F90"
    "${top_srcdir}/radfft.f"
    "${top_srcdir}/radial.f"
    "${top_srcdir}/reclat.f"
    "${top_srcdir}/read_xc_info.F"
    "${top_srcdir}/siesta_geom.F90"
    "${top_srcdir}/sorting.f"
    "${top_srcdir}/xml.f"
     ${top_srcdir}/m_spin.F90
     ${top_srcdir}/m_vee_integrals.F90

)

if(NOT LIBGRIDXC_USES_PROCEDURE_POINTER)
  #
  # We might need to resolve the 'timer' reference in gridxc_wrappers
  #
  list(APPEND sources  local_timer.f90)
endif()

siesta_add_library(${PROJECT_NAME}.aux-lib-basis OBJECT ${sources})

#
# lapack is needed for m_filter...
target_link_libraries(
  ${PROJECT_NAME}.aux-lib-basis
  PRIVATE
  libfdf::libfdf
  libgridxc::libgridxc
  libpsml::libpsml
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libpsop
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
  LAPACK::LAPACK
 )


siesta_add_executable(${PROJECT_NAME}.gen-basis gen-basis.F)
target_link_libraries(
  ${PROJECT_NAME}.gen-basis
  PRIVATE
  ${PROJECT_NAME}.aux-lib-basis
  libfdf::libfdf
  libgridxc::libgridxc
  libpsml::libpsml
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libpsop
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  LAPACK::LAPACK
 )


siesta_add_executable(${PROJECT_NAME}.ioncat ioncat.f
  ${top_srcdir}/atmfuncs.f
  ${top_srcdir}/spher_harm.f
  ${top_srcdir}/m_getopts.f90
  ${top_srcdir}/chkdim.f
)

target_link_libraries(
  ${PROJECT_NAME}.ioncat
  PRIVATE
  ${PROJECT_NAME}.aux-lib-basis
  libfdf::libfdf
  libgridxc::libgridxc
  libpsml::libpsml
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libpsop
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  LAPACK::LAPACK
  ${PROJECT_NAME}.mpi_macros
 )



if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.gen-basis ${PROJECT_NAME}.ioncat
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

  #
  # Copy also ionplot.sh (with executable bits?)
  #
  install(
    FILES
    "ionplot.sh"
    PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
                GROUP_EXECUTE GROUP_READ
                WORLD_EXECUTE WORLD_READ
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
