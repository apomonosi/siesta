#
# There are two targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
# Since the top objects are compiled separately, their 'sys' module
# must be the top one. Here only external 'die' et al are bundled in
# the 'sys_externals.f90' file.
#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

siesta_add_library(${PROJECT_NAME}.vibra_top_objs OBJECT
   ${top_src_dir}/io.f
   ${top_src_dir}/m_io.f
   ${top_src_dir}/precision.F
   ${top_src_dir}/reclat.f
)

target_link_libraries(
    ${PROJECT_NAME}.vibra_top_objs
    PRIVATE
    ${PROJECT_NAME}.libsys
    )

siesta_add_executable(${PROJECT_NAME}.fcbuild
    fcbuild.f
    recoor.f
)

siesta_add_executable(${PROJECT_NAME}.vibra
   hermdp.F
   klines.f
   outbands.f
   recoor.f
   vibra.F
)

target_link_libraries(${PROJECT_NAME}.fcbuild
  PRIVATE
    ${PROJECT_NAME}.vibra_top_objs
    libfdf::libfdf
    ${PROJECT_NAME}.libsys
    ${PROJECT_NAME}.libunits
    )
target_link_libraries(${PROJECT_NAME}.vibra
  PRIVATE
    ${PROJECT_NAME}.vibra_top_objs
    libfdf::libfdf
    ${PROJECT_NAME}.libsys
    ${PROJECT_NAME}.libunits
    LAPACK::LAPACK
    )
    
if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.fcbuild ${PROJECT_NAME}.vibra
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
