# This is probably the simplest solution
# to add the C++ runtime to the linking statement
# Useful for, e.g., PEXSI
#
siesta_add_library(${PROJECT_NAME}.cxx_dummy_lib
  NO_LINKER_FLAGS
  OBJECT
    dummy.cpp
  )

target_include_directories(${PROJECT_NAME}.cxx_dummy_lib
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}"
)
