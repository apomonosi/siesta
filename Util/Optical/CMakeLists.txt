siesta_add_executable( ${PROJECT_NAME}.optical optical_cond.f )
siesta_add_executable( ${PROJECT_NAME}.optical_input optical_input.f )

target_link_libraries(${PROJECT_NAME}.optical
  PRIVATE ${PROJECT_NAME}.libunits ${PROJECT_NAME}.libsys)
target_link_libraries(${PROJECT_NAME}.optical_input
  PRIVATE ${PROJECT_NAME}.libunits ${PROJECT_NAME}.libsys)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.optical ${PROJECT_NAME}.optical_input
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

