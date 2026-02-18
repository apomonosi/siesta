siesta_add_executable(${PROJECT_NAME}.pdosxml
     m_orbital_chooser.f90 m_pdos.f90 pdosxml.f90)

target_link_libraries(${PROJECT_NAME}.pdosxml PRIVATE xmlf90::xmlf90)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.pdosxml
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

