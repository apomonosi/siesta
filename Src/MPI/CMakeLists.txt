#
# When to use the legacy interfaces:
#
#    When mpi_f08 is not available or wanted and
#    the compiler is able to process them. If it is not
#    (e.g. some versions of nvhpc), SIESTA_WITH_MPI_EMPTY_INTERFACES
#    should be defined.
#
#

# Create the MPI mpi_macros.f interface library.
# This just simplifies the usage of the
#  #include "mpi_macros.f"
# in the files.
add_library(${PROJECT_NAME}.mpi_macros INTERFACE)

if( SIESTA_WITH_MPI )
  # Real mpi_macros.f
  target_include_directories(${PROJECT_NAME}.mpi_macros
    INTERFACE
    ${PROJECT_SOURCE_DIR}/Src/MPI
  )
else()
  # fake mpi_macros.f
  target_include_directories(${PROJECT_NAME}.mpi_macros
    INTERFACE
    ${PROJECT_SOURCE_DIR}/Src/MPI/stub
  )

  # Quick-return to not do anything
  return()
endif()


if(SIESTA_WITH_MPI_INTERFACES STREQUAL "legacy")
    set(SIESTA_USES_MPI_LEGACY_INTERFACES ON)
else()
    set(SIESTA_USES_MPI_LEGACY_INTERFACES OFF)
endif()

message(STATUS "SIESTA_USES_MPI_LEGACY_INTERFACES: ${SIESTA_USES_MPI_LEGACY_INTERFACES}")

if (NOT SIESTA_USES_MPI_LEGACY_INTERFACES)

  siesta_add_library(${PROJECT_NAME}.mpi
    NO_LINKER_FLAGS
    OBJECT
    mpi_siesta.F90
  )

else()

  siesta_add_library(${PROJECT_NAME}.mpi
    NO_LINKER_FLAGS
    OBJECT

    # Sources
    mpi_siesta.F90
    mpi_interfaces.F90
    mpi__include.f90
    Interfaces_integer.f90
    Interfaces.f90
  )
  # We need to find "Interfaces.f90", which is generated in the build
  #("binary") directory, so we need to add this directory to the list

  target_include_directories(${PROJECT_NAME}.mpi
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

endif()

string(TOUPPER ${SIESTA_WITH_MPI_INTERFACES} INTERFACE_upper)
set(VAR_MPI_INTERFACE "MPI_INTERFACE_${INTERFACE_upper}")
message(STATUS "Variable for preprocessing MPI interfaces: ${VAR_MPI_INTERFACE}")

target_compile_definitions(
  ${PROJECT_NAME}.mpi
  PUBLIC
  "$<$<BOOL:${SIESTA_WITH_GRID_SP}>:GRID_SP>"
  "${VAR_MPI_INTERFACE}"
)

target_link_libraries(${PROJECT_NAME}.mpi
  PUBLIC MPI::MPI_Fortran
  INTERFACE ${PROJECT_NAME}.mpi_macros
)

# So that clients of ${PROJECT_NAME}.mpi can find its modules:
target_include_directories(${PROJECT_NAME}.mpi
  INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR} # for created files
)


if (SIESTA_USES_MPI_LEGACY_INTERFACES)
  # Handling of interface-generation script
  # These generation scripts do not work with the suffix for the executables
  siesta_suffix(
    NEW ""
    VARIABLES EXECUTABLE)


  # Create integer kinds and interfaces
  siesta_add_executable(${PROJECT_NAME}.mpi_int_explorer
    # Do not default to build this one, it is only
    # a temporary build dependency
    EXCLUDE_FROM_ALL
    int_explorer.f90
  )


  list(JOIN SIESTA_INTEGER_KINDS " " SIESTA_INTEGER_KINDS_STR)
  add_custom_command(
    OUTPUT Interfaces_integer.f90
    # Dependency targets
    DEPENDS ${PROJECT_NAME}.mpi_int_explorer
    COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/generate_integer.sh" ${SIESTA_INTEGER_KINDS_STR}
    VERBATIM
    COMMENT "Creating MPI interfaces and modules for integer kinds"
  )


  # Create real kinds and interfaces
  siesta_add_executable(${PROJECT_NAME}.mpi_kind_explorer
    EXCLUDE_FROM_ALL
    kind_explorer.f90
  )

  list(JOIN SIESTA_REAL_KINDS " " SIESTA_REAL_KINDS_STR)
  add_custom_command(
    OUTPUT Interfaces.f90
    DEPENDS ${PROJECT_NAME}.mpi_kind_explorer
    COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/generate_real.sh" ${SIESTA_REAL_KINDS_STR}
    VERBATIM
    COMMENT "Creating MPI interfaces and modules for real, character and logical"
  )

  # Restore the old suffix before using this one.
  siesta_suffix(
    POP
    VARIABLES EXECUTABLE
  )
endif()


