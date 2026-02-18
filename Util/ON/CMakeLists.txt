if (SIESTA_WITH_NETCDF)

  siesta_add_executable(${PROJECT_NAME}.lwf2cdf lwf2cdf.F90 )

  if( SIESTA_INSTALL )
    install(
      TARGETS ${PROJECT_NAME}.lwf2cdf
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
  endif()

endif()


