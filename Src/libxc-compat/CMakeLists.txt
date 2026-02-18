siesta_add_library(${PROJECT_NAME}.libxc-trans
  NO_LINKER_FLAGS
  OBJECT
    m_libxc_compat.f90
    m_libxc_list.F90
    m_libxc_sxc_translation.F90
    m_siestaxc_list.F90
)

target_include_directories(
  ${PROJECT_NAME}.libxc-trans
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

