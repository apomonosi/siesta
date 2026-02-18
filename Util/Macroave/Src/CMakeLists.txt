set(top_srcdir "${PROJECT_SOURCE_DIR}/Src" )

set(sources

  defs_basis.f90
  defs_common.f90
  io.f
  iorho.f
  macroave.f
  thetaft.f
  surpla.f
  volcel.f
  hdr_io.f90
  dpnint1.f90

)

list(
  APPEND
  sources

  ${top_srcdir}/interpolation.f90
  ${top_srcdir}/m_fft_gpfa.F
  ${top_srcdir}/precision.F

)

siesta_add_executable(${PROJECT_NAME}.macroave ${sources})
siesta_add_executable(${PROJECT_NAME}.permute permute.F)

target_link_libraries(${PROJECT_NAME}.macroave
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

if( SIESTA_INSTALL )
  install(
    TARGETS
    ${PROJECT_NAME}.macroave
    ${PROJECT_NAME}.permute
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

