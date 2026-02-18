set(qmmmtop_srcdir "${PROJECT_SOURCE_DIR}/Src")

set(libsqmmm_sources
  siestaqmmm.F90
  siestaqmmm_m.F90
)

list(
  APPEND
  libsqmmm_sources
  ${qmmmtop_srcdir}/fsockets.f90
  ${qmmmtop_srcdir}/precision.F
  ${qmmmtop_srcdir}/sockets.c
)

add_library(${PROJECT_NAME}-libsiestaqmmm "${libsqmmm_sources}")
set_target_properties(
  ${PROJECT_NAME}-libsiestaqmmm
  PROPERTIES
  VERSION "${PROJECT_VERSION}"
  OUTPUT_NAME "siestaqmmm"
)

install(
  TARGETS ${PROJECT_NAME}-libsiestaqmmm
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

target_compile_definitions(
    ${PROJECT_NAME}-libsiestaqmmm
    PUBLIC
    "$<$<BOOL:${SIESTA_WITH_GRID_SP}>:GRID_SP>"
)

target_link_libraries(
  ${PROJECT_NAME}-libsiestaqmmm
  PRIVATE
  $<$<BOOL:${SIESTA_WITH_OPENMP}>:OpenMP::OpenMP_Fortran>
  libfdf::libfdf
)
