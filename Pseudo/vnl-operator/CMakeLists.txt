set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

add_library(${PROJECT_NAME}.psop_top_lib OBJECT
    "${top_srcdir}/m_getopts.f90"
    "${top_srcdir}/m_uuid.f90"
    "${top_srcdir}/periodic_table.f"
)

set(psop_sources
  check_grid.f90
  dpnint.f90
  m_kb.f90
  psop.f90
  psop_options.f90
)

if(NOT LIBGRIDXC_USES_PROCEDURE_POINTER)
  #
  # We might need to resolve the 'timer' reference in gridxc_wrappers
  #
  list(APPEND psop_sources  local_timer.f90)
endif()

siesta_add_executable(
  ${PROJECT_NAME}.psop
  ${psop_sources}
)

# For semicore_info_froyen...
target_link_libraries(${PROJECT_NAME}.psop_top_lib
  PRIVATE
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libsys
 )
 
target_link_libraries(
  ${PROJECT_NAME}.psop
  PRIVATE
  ${PROJECT_NAME}.psop_top_lib
  libgridxc::libgridxc
  libpsml::libpsml
  xmlf90::xmlf90
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libpsop
  ${PROJECT_NAME}.libxc-trans
  ${PROJECT_NAME}.libsys
 )

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.psop
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
