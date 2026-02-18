#
siesta_add_executable(${PROJECT_NAME}.sies2arc
   sies2arc.f linepro.f wtof.f cell.f uncell.f
)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.sies2arc
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

