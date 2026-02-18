set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

set(sources
    ts2ts.f90 )

list(APPEND sources

   ${top_srcdir}/precision.F
   ${top_srcdir}/parallel.F
   ${top_srcdir}/alloc.F90
   ${top_srcdir}/m_io.f
)

siesta_add_executable(${PROJECT_NAME}.ts2ts ${sources})
target_link_libraries(${PROJECT_NAME}.ts2ts
  PRIVATE
  ${PROJECT_NAME}.libsys
  ${PROJECT_NAME}.libunits
  ${PROJECT_NAME}.mpi_macros
	libfdf::libfdf
)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.ts2ts
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
