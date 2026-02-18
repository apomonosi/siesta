siesta_add_executable(${PROJECT_NAME}.grid1d
  grid1d.f )
siesta_add_executable(${PROJECT_NAME}.grid2d
  grid2d.f )

target_link_libraries(${PROJECT_NAME}.grid1d
  PRIVATE
  ${PROJECT_NAME}.libunits
  )
target_link_libraries(${PROJECT_NAME}.grid2d
  PRIVATE
  ${PROJECT_NAME}.libunits
  )


if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.grid1d ${PROJECT_NAME}.grid2d
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
