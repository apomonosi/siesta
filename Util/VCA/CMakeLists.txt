#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

siesta_add_library(${PROJECT_NAME}.vca_objs OBJECT
   ${top_src_dir}/interpolation.f90
   ${top_src_dir}/m_io.f
   ${top_src_dir}/io.f
   ${top_src_dir}/periodic_table.f
   ${top_src_dir}/precision.F
)
target_link_libraries(${PROJECT_NAME}.vca_objs
  PRIVATE
  ${PROJECT_NAME}.libsys)

siesta_add_executable( ${PROJECT_NAME}.mixps mixps.f )
siesta_add_executable( ${PROJECT_NAME}.fractional fractional.f)


target_link_libraries(${PROJECT_NAME}.mixps
  PRIVATE
  ${PROJECT_NAME}.vca_objs
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libsys
  )
target_link_libraries(${PROJECT_NAME}.fractional
  PRIVATE
  ${PROJECT_NAME}.vca_objs
  ${PROJECT_NAME}.libncps
  ${PROJECT_NAME}.libsys
  )
    
if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.mixps ${PROJECT_NAME}.fractional
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
