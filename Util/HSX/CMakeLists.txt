#
# Go the extra mile and treat this as a library
# that users can link to in their own programs.
#
siesta_add_library(${PROJECT_NAME}.libhsx
  OUTPUT_NAME hsx
  hsx_m.f90)

set_target_properties(${PROJECT_NAME}.libhsx
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)

target_include_directories(
  ${PROJECT_NAME}.libhsx
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

#
# Now the utilities
# 
siesta_add_executable(${PROJECT_NAME}.hsx2hs hsx2hs.f90)
siesta_add_executable(${PROJECT_NAME}.hs2hsx hs2hsx.f90)

target_link_libraries(${PROJECT_NAME}.hsx2hs PRIVATE ${PROJECT_NAME}.libhsx)
target_link_libraries(${PROJECT_NAME}.hs2hsx PRIVATE ${PROJECT_NAME}.libhsx)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.libhsx
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  install(
    DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
  install(
    TARGETS
      ${PROJECT_NAME}.hsx2hs ${PROJECT_NAME}.hs2hsx
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
