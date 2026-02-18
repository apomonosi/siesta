set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")
set(grid_srcdir "${PROJECT_SOURCE_DIR}/Util/Grid")

# To avoid "multiple rules" in Ninja
# Maybe turn into global-utility targets

siesta_add_library(${PROJECT_NAME}.stm_gridfunc   ${grid_srcdir}/m_gridfunc.F90)

siesta_add_executable( ${PROJECT_NAME}.plstm 
  plstm.f90
  ${top_srcdir}/m_getopts.f90
)

siesta_add_executable( ${PROJECT_NAME}.plsts
  plsts.f
)

target_link_libraries(${PROJECT_NAME}.plstm PRIVATE ${PROJECT_NAME}.stm_gridfunc)
target_link_libraries(${PROJECT_NAME}.plsts PRIVATE ${PROJECT_NAME}.stm_gridfunc)

if( SIESTA_INSTALL )
  install(
    TARGETS ${PROJECT_NAME}.plstm ${PROJECT_NAME}.plsts
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

