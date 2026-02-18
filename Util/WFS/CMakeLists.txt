#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

#-----------------------------------------
siesta_add_executable(${PROJECT_NAME}.readwfx
   ${top_src_dir}/m_getopts.f90
   readwfx.f90
)

siesta_add_executable(${PROJECT_NAME}.readwf readwf.f)  
siesta_add_executable(${PROJECT_NAME}.wfs2wfsx  wfs2wfsx.f)
siesta_add_executable(${PROJECT_NAME}.wfsx2wfs  wfsx2wfs.f)

if( SIESTA_INSTALL )
  install(
    TARGETS
      ${PROJECT_NAME}.readwf ${PROJECT_NAME}.readwfx
      ${PROJECT_NAME}.wfs2wfsx ${PROJECT_NAME}.wfsx2wfs
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif()

if (SIESTA_WITH_NETCDF)

  siesta_add_executable(${PROJECT_NAME}.wfsnc2wfsx wfsnc2wfsx.F90)

  target_link_libraries(${PROJECT_NAME}.wfsnc2wfsx
    PRIVATE
    NetCDF::NetCDF_Fortran
    ${PROJECT_NAME}.libunits
    )
  target_compile_definitions(${PROJECT_NAME}.wfsnc2wfsx  PRIVATE CDF )
		     
  if( SIESTA_INSTALL )
    install(
      TARGETS ${PROJECT_NAME}.wfsnc2wfsx
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
  endif()

endif(SIESTA_WITH_NETCDF)
