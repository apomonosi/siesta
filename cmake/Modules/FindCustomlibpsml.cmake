include(SiestaFindPackage)

Siesta_find_package(libpsml
  REQUIRED
  GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/libpsml.git"
  GIT_TAG "master"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/libpsml
  )

# We might know about this feature from an explicit record in the source

if( "${libpsml_FOUND_METHOD}" STREQUAL "cmake" OR
    "${libpsml_FOUND_METHOD}" STREQUAL "pkgconf")

  include(CheckFortranSourceCompiles)

  # Figure out whether psml uses procedure pointers, or not
  set(CMAKE_REQUIRED_LIBRARIES libpsml::libpsml)
  check_fortran_source_compiles("use m_psml, only: ps_set_error_handler; end"
    LIBPSML_HAS_ERROR_PROCEDURE_POINTER SRC_EXT F90)
  unset(CMAKE_REQUIRED_LIBRARIES)

else()

  if(LIBPSML_HAS_ERROR_PROCEDURE_POINTER)
     message(STATUS "  libpsml (submodule/fetch) has the 'procedure pointer' feature")
  endif()

endif()



if( LIBPSML_HAS_ERROR_PROCEDURE_POINTER )
  set(LIBPSML_USES_PROCEDURE_POINTER TRUE CACHE BOOL "Whether libpsml uses procedure pointers or not")
else()
  set(LIBPSML_USES_PROCEDURE_POINTER FALSE CACHE BOOL "Whether libpsml uses procedure pointers or not")
endif()
mark_as_advanced(LIBPSML_USES_PROCEDURE_POINTER)
