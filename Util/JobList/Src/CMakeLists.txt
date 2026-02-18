# To avoid "multiple rules" in Ninja
# Maybe turn into global-utility targets

siesta_add_library(${PROJECT_NAME}.aux_joblist
  jobList.f90 posix_calls.f90)

siesta_add_executable(${PROJECT_NAME}.countJobs  countJobs.f90 )
siesta_add_executable(${PROJECT_NAME}.runJobs    runJobs.f90 )
siesta_add_executable(${PROJECT_NAME}.getResults getResults.f90 )
siesta_add_executable(${PROJECT_NAME}.horizontal horizontal.f90 )

target_link_libraries(${PROJECT_NAME}.countJobs PRIVATE ${PROJECT_NAME}.aux_joblist)
target_link_libraries(${PROJECT_NAME}.runJobs PRIVATE ${PROJECT_NAME}.aux_joblist)
target_link_libraries(${PROJECT_NAME}.getResults PRIVATE ${PROJECT_NAME}.aux_joblist)
target_link_libraries(${PROJECT_NAME}.horizontal PRIVATE ${PROJECT_NAME}.aux_joblist)
		     
if( SIESTA_INSTALL )
  install(TARGETS
    ${PROJECT_NAME}.countJobs
    ${PROJECT_NAME}.runJobs
    ${PROJECT_NAME}.getResults
    ${PROJECT_NAME}.horizontal
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} 
  )
endif()
