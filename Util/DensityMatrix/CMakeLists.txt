add_library(${PROJECT_NAME}.objs_getopts OBJECT
            ${PROJECT_SOURCE_DIR}/Src/m_getopts.f90 )
	    
siesta_add_executable(${PROJECT_NAME}.dmfilter dmfilter.f90)
siesta_add_executable(${PROJECT_NAME}.dmbs2dm dmbs2dm.F90 )
target_link_libraries(${PROJECT_NAME}.dmbs2dm PRIVATE ${PROJECT_NAME}.objs_getopts)
siesta_add_executable(${PROJECT_NAME}.dmUnblock dmUnblock.F90)
target_link_libraries(${PROJECT_NAME}.dmUnblock PRIVATE ${PROJECT_NAME}.objs_getopts)
siesta_add_executable(${PROJECT_NAME}.dm_noncol_sign_flip4 dm_noncol_sign_flip4.f90)

if( SIESTA_INSTALL )
  install(
    TARGETS
      ${PROJECT_NAME}.dmfilter
      ${PROJECT_NAME}.dmbs2dm
      ${PROJECT_NAME}.dmUnblock
      ${PROJECT_NAME}.dm_noncol_sign_flip4
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

if(SIESTA_WITH_NETCDF)

  siesta_add_executable(${PROJECT_NAME}.dm2cdf dm2cdf.F90)
  siesta_add_executable(${PROJECT_NAME}.cdf2dm cdf2dm.F90)

  if( SIESTA_INSTALL )
    install(
      TARGETS
        ${PROJECT_NAME}.dm2cdf
        ${PROJECT_NAME}.cdf2dm
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
  endif()
  
endif(SIESTA_WITH_NETCDF)

