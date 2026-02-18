#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

add_library(${PROJECT_NAME}.bands_top_objs OBJECT
   ${top_src_dir}/m_getopts.f90
)

siesta_add_executable(${PROJECT_NAME}.gnubands
    gnubands.f90
)

siesta_add_executable(${PROJECT_NAME}.eigfat2plot
   eigfat2plot.f90
)

target_link_libraries(${PROJECT_NAME}.gnubands
  PRIVATE
  ${PROJECT_NAME}.bands_top_objs
  )
target_link_libraries(${PROJECT_NAME}.eigfat2plot
  PRIVATE
  ${PROJECT_NAME}.bands_top_objs
  )

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.gnubands ${PROJECT_NAME}.eigfat2plot
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
