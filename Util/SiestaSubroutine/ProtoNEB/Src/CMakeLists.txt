if( NOT (SIESTA_WITH_MPI AND NOT SIESTA_WITH_NO_MPI_INTERFACES) )
  # this package should not be added
  message(WARNING "${CMAKE_CURRENT_SOURCE_DIR} skipped due to missing MPI support")
  return()
endif()

set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

# ProtoNEB: Uses the MPI dispatch
#
siesta_add_executable(${PROJECT_NAME}.protoNEB
  NAMESPACE_TARGET protoNEB

   ${top_srcdir}/reinit_m.F90
   ${top_srcdir}/siesta_cmlsubs.F90
   ${top_srcdir}/siesta_init.F
   ${top_srcdir}/m_io_yaml.F90
   ${top_srcdir}/siesta_end.F
   ${top_srcdir}/fsiesta_mpi.F90
   ${top_srcdir}/init_output.f90

    protoNEB.F90
)

target_link_libraries(
  ${PROJECT_NAME}.protoNEB
  PRIVATE
  ${PROJECT_NAME}.libsiesta
  )

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.protoNEB
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

