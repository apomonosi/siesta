siesta_add_executable(${PROJECT_NAME}.orbmol_proj orbmol_proj.f90)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.orbmol_proj
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

