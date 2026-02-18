#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

siesta_add_executable(${PROJECT_NAME}.psml2psf
   ${top_src_dir}/m_getopts.f90

   psml2psf.f90
)

target_link_libraries(
  ${PROJECT_NAME}.psml2psf
  PRIVATE
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libsys
  libpsml::libpsml
  )

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.psml2psf
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
