siesta_add_library(
  ${PROJECT_NAME}.libsys
  NO_LINKER_FLAGS
  OBJECT
    sys.F90
    external_entries.F90
  )

target_include_directories(
  ${PROJECT_NAME}.libsys
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}
)

