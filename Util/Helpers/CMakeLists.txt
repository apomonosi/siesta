#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

siesta_add_executable(${PROJECT_NAME}.get_chem_labels
   ${top_src_dir}/parallel.F
   ${top_src_dir}/chemical.f
   ${top_src_dir}/atom_options.F90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/io.f
   ${top_src_dir}/alloc.F90
   ${top_src_dir}/memory.F
   ${top_src_dir}/memory_log.F90
   ${top_src_dir}/periodic_table.f
   ${top_src_dir}/precision.F
   ${top_src_dir}/pxf.F90
   get_chem_labels.f90
   local_timer.f90
)

target_link_libraries(${PROJECT_NAME}.get_chem_labels
  PRIVATE
	libfdf::libfdf
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
  )


install(
  TARGETS ${PROJECT_NAME}.get_chem_labels
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# This executable will not be installed.
