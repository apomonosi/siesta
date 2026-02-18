#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

siesta_add_library(${PROJECT_NAME}.sockets_objs OBJECT
   ${top_src_dir}/fsockets.f90
   ${top_src_dir}/sockets.c
)

siesta_add_executable( ${PROJECT_NAME}.f2fmaster f2fmaster.f90 )
siesta_add_executable( ${PROJECT_NAME}.f2fslave  f2fslave.f90 )

target_link_libraries(${PROJECT_NAME}.f2fmaster PRIVATE ${PROJECT_NAME}.sockets_objs )
target_link_libraries(${PROJECT_NAME}.f2fslave PRIVATE ${PROJECT_NAME}.sockets_objs )
    
if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.f2fmaster ${PROJECT_NAME}.f2fslave
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

