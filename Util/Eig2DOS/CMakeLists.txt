set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")
siesta_add_executable(${PROJECT_NAME}.Eig2DOS
   Eig2DOS.f90
   ${top_src_dir}/m_getopts.f90
)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.Eig2DOS
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
