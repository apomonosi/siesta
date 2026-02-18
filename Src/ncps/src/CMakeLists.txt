siesta_add_library(
  ${PROJECT_NAME}.libncps
  NO_LINKER_FLAGS
  STATIC

    interpolation.f90
    m_ncps.f90
    m_ncps_froyen_ps_t.f90
    m_ncps_froyen_reader.f
    m_ncps_psml_plugins.f90
    m_ncps_reader.f
    m_ncps_translators.f90
    m_ncps_utils.f90
    m_ncps_writers.f90
    env_utils_m.f90
    search_ps_m.F90
  )

target_link_libraries(
  ${PROJECT_NAME}.libncps
  PUBLIC
    ${PROJECT_NAME}.libxc-trans
    ${PROJECT_NAME}.libsys
    libpsml::libpsml
)
target_include_directories(
  ${PROJECT_NAME}.libncps
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

