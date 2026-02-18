#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

siesta_add_executable(${PROJECT_NAME}.fdf2grimme
   ${top_src_dir}/parallel.F
   ${top_src_dir}/chemical.f
   ${top_src_dir}/atom_options.F90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/io.f
   ${top_src_dir}/periodic_table.f
   ${top_src_dir}/precision.F
   grimme.f90
   local_timer.f90
)

target_link_libraries(${PROJECT_NAME}.fdf2grimme
  PRIVATE
  libfdf::libfdf
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.fdf2grimme
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
