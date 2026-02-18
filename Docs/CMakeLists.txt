#
# Logic taken from SIRIUS distribution
#
#
# Doxygen version of the documentation
#
if( NOT SIESTA_BUILD_DOCS )
  return()
endif()

find_package(Doxygen)

if( DOXYGEN_FOUND )

  set(DOXYGEN_CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/doxygen.cfg")
  # configure doxygen.cfg.in and adapt paths
  configure_file(doxygen.cfg.in ${DOXYGEN_CONFIG_FILE} @ONLY IMMEDIATE)

  set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  # Add a target and build doc with all targets
  add_custom_target(doc-doxygen ALL
     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG_FILE}
     WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIR}
   )
  install(
    DIRECTORY ${DOXYGEN_OUTPUT_DIR}/html/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/doc/doxygen
  )
endif()

#
# FORD version of the documentation
#
# Some modern versions of FORD are unable to process correctly the
# inter-page links. Treat this as work-in-progress
#
find_program(FORD_EXE ford
  DOC "path to the ford executable (required to generate the documentation)"
)
if( FORD_EXE )

  # Copy the FORD project-file into the build directory
  set(FORD_PROJECT_FILE "${CMAKE_CURRENT_BINARY_DIR}/ford-siesta.md")
  set(FORD_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc-ford-build")
  configure_file(ford.cfg.in "${FORD_PROJECT_FILE}" @ONLY)

  add_custom_target(doc-ford
    COMMENT "Generating API documentation"
    COMMAND "${FORD_EXE}" "${FORD_PROJECT_FILE}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    VERBATIM
  )
  install(
    DIRECTORY "${FORD_OUTPUT_DIR}/"
    DESTINATION ${CMAKE_INSTALL_PREFIX}/doc/ford
  )

endif()
