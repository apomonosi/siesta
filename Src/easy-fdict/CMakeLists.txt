# This internal library uses a "pre-processed" file from the fdict
# distribution, containing the needed variable types.
  
siesta_add_library(
  ${PROJECT_NAME}.libfdict
  NO_LINKER_FLAGS
  OBJECT
    variable.F90
    dictionary.f90
  )

target_include_directories(
  ${PROJECT_NAME}.libfdict
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}"
)

