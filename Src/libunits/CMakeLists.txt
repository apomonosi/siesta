siesta_add_library(
  ${PROJECT_NAME}.libunits
  NO_LINKER_FLAGS
  STATIC
    units.F90
    units_legacy.f90
    units_codata2018.f90
    units_common.f90
  )

target_link_libraries(
  ${PROJECT_NAME}.libunits
  PUBLIC
  ${PROJECT_NAME}.libsys
)


# This is duplicated in ${PROJECT_SOURCE_DIR}/CMakeLists.txt
# The reason is that we might remove the legacy code in all codes
# and fully resort to only handling this in units* files in this library.
# That would be the longterm goal.

message(STATUS "Value of SIESTA_WITH_UNIT_CONVENTION in Src/libunits: ${SIESTA_WITH_UNIT_CONVENTION}")

# Lower case the units
string(TOLOWER "${SIESTA_WITH_UNIT_CONVENTION}" siesta_unit_convention)
if("${siesta_unit_convention}" STREQUAL "original" OR
   "${siesta_unit_convention}" STREQUAL "legacy")
  target_compile_definitions(${PROJECT_NAME}.libunits
    PUBLIC
    SIESTA__UNITS_ORIGINAL)
elseif("${siesta_unit_convention}" STREQUAL "codata2018")
  # pass - do nothing - the default
else()
  message(FATAL_ERROR "Unknown argument for SIESTA_WITH_UNIT_CONVENTION; not one of [original|legacy,codata2018]")
endif()
message(STATUS "Siesta using units: ${siesta_unit_convention}")

# Ensure module files are added to the inclusion
target_include_directories(
  ${PROJECT_NAME}.libunits
  INTERFACE "${CMAKE_CURRENT_BINARY_DIR}")

