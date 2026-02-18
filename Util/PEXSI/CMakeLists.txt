if( NOT SIESTA_WITH_PEXSI )
  return()
endif()

siesta_add_executable(${PROJECT_NAME}.intdos2eig intdos2eig.f90)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.intdos2eig
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

