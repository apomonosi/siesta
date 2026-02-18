siesta_add_library(
  ${PROJECT_NAME}.libwxml
  NO_LINKER_FLAGS
  OBJECT
     
    flib_wcml.f90
    flib_wstml.f90
    flib_wxml.f90
    m_wcml_coma.f90
    m_wxml_array_str.f90
    m_wxml_buffer.f90
    m_wxml_core.f90
    m_wxml_dictionary.f90
    m_wxml_elstack.f90
    m_wxml_error.f90
    m_wxml_escape.f90
    m_wxml_overloads.f90
    m_wxml_text.F90
  )

target_include_directories(
  ${PROJECT_NAME}.libwxml
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

