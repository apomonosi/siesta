siesta_add_executable(
  ${PROJECT_NAME}.test-ncps

   die.f90
   m_getopts.f90
   psml_die.F90
   test_ncps.f90
)

### install(
###   TARGETS test-ncps
###   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
###   )

target_link_libraries(${PROJECT_NAME}.test-ncps
   ${PROJECT_NAME}.libncps
   libpsml::libpsml
   xmlf90::xmlf90
)

