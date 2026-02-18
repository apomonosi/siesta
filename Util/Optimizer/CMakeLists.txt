siesta_add_library(${PROJECT_NAME}.aux_optim
  OBJECT
    io.f
    parse.f
    minimizer.f90
    vars_module.f90
    precision.f90
    system.f90
    sys.f
)    

siesta_add_executable(${PROJECT_NAME}.simplex
    amoeba.f
    simplex.f90
)

siesta_add_executable(${PROJECT_NAME}.swarm
    swarm.f90
)

target_link_libraries(${PROJECT_NAME}.aux_optim PUBLIC
  $<$<BOOL:${SIESTA_WITH_OPENMP}>:OpenMP::OpenMP_Fortran>
  )

target_link_libraries(${PROJECT_NAME}.swarm PRIVATE ${PROJECT_NAME}.aux_optim)
target_link_libraries(${PROJECT_NAME}.simplex PRIVATE ${PROJECT_NAME}.aux_optim)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.simplex ${PROJECT_NAME}.swarm
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

