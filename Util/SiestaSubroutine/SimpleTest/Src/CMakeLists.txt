#
# FIXME: The names in this directory are very confusing
#        Only a subset of the executables are built
#
set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")


#
# Create auxiliary libraries, with their own module directories,
# to satisfy Ninja. This is also good practice in general
#
# -- Sockets
siesta_add_library(${PROJECT_NAME}.sockets-dispatch
  NO_LINKER_FLAGS
  OBJECT
    ${top_srcdir}/fsiesta_sockets.F90
    ${top_srcdir}/fsockets.f90
    ${top_srcdir}/sockets.c
    ${top_srcdir}/posix_calls.f90
)
set_target_properties(${PROJECT_NAME}.sockets-dispatch
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sockets-include")

target_include_directories(${PROJECT_NAME}.sockets-dispatch
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/sockets-include")

# -- Pipes
siesta_add_library(${PROJECT_NAME}.pipes-dispatch
  NO_LINKER_FLAGS
  OBJECT
    ${top_srcdir}/fsiesta_pipes.F90
    ${top_srcdir}/pxf.F90
    ${top_srcdir}/posix_calls.f90
)
set_target_properties(${PROJECT_NAME}.pipes-dispatch
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pipes-include")
target_include_directories(${PROJECT_NAME}.pipes-dispatch
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/pipes-include")

if (SIESTA_WITH_MPI)

# -- MPI dispatch -- uses 'libsiesta'
siesta_add_library(${PROJECT_NAME}.mpi-dispatch
  NO_LINKER_FLAGS
  OBJECT
    ${top_srcdir}/fsiesta_mpi.F90
)
target_link_libraries(${PROJECT_NAME}.mpi-dispatch PUBLIC ${PROJECT_NAME}.libsiesta)

set_target_properties(${PROJECT_NAME}.mpi-dispatch
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mpi-include")
target_include_directories(${PROJECT_NAME}.mpi-dispatch
  INTERFACE
  "${CMAKE_CURRENT_BINARY_DIR}/mpi-include")

# ----
#
# Phonons: Uses the MPI dispatch
#
siesta_add_executable(${PROJECT_NAME}.phonons
  phonons.F90
)
target_link_libraries(${PROJECT_NAME}.phonons PRIVATE ${PROJECT_NAME}.mpi-dispatch)

if(MPI_Fortran_HAVE_F08_MODULE)
  siesta_add_executable(${PROJECT_NAME}.phonons-f08
  phonons.F90
  )
  target_link_libraries(${PROJECT_NAME}.phonons-f08 PRIVATE ${PROJECT_NAME}.mpi-dispatch)
  target_compile_definitions(${PROJECT_NAME}.phonons-f08 PRIVATE CLIENT_USES_F08)
endif()

#
# Test in which not all ranks are involved in siesta_launch
#
siesta_add_executable(${PROJECT_NAME}.comm_test
  comm_test.f90)
target_link_libraries(${PROJECT_NAME}.comm_test PRIVATE ${PROJECT_NAME}.mpi-dispatch)

# Simple MPI-aware driver
siesta_add_executable(${PROJECT_NAME}.mpi_driver
  simple_mpi_parallel.f90
)
target_link_libraries(${PROJECT_NAME}.mpi_driver PRIVATE ${PROJECT_NAME}.mpi-dispatch)


# Sockets interface to parallel-version of Siesta
siesta_add_executable(${PROJECT_NAME}.sockets_parallel
  simple_parallel.f90)
target_link_libraries(${PROJECT_NAME}.sockets_parallel PRIVATE ${PROJECT_NAME}.sockets-dispatch)


# Pipes interface to parallel-version of Siesta
siesta_add_executable(${PROJECT_NAME}.pipes_parallel
  simple_parallel.f90
)
target_link_libraries(${PROJECT_NAME}.pipes_parallel PRIVATE ${PROJECT_NAME}.pipes-dispatch)


if( SIESTA_INSTALL )
  install(
    TARGETS
    ${PROJECT_NAME}.phonons ${PROJECT_NAME}.pipes_parallel
      ${PROJECT_NAME}.sockets_parallel ${PROJECT_NAME}.mpi_driver
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
  if(MPI_Fortran_HAVE_F08_MODULE)
    install(
      TARGETS
      ${PROJECT_NAME}.phonons-f08
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      )
  endif()
endif()
endif(SIESTA_WITH_MPI)

siesta_add_executable(${PROJECT_NAME}.sockets_serial
  simple_serial.f90)
target_link_libraries(${PROJECT_NAME}.sockets_serial PRIVATE ${PROJECT_NAME}.sockets-dispatch)

siesta_add_executable(${PROJECT_NAME}.pipes_serial
  simple_serial.f90)
target_link_libraries(${PROJECT_NAME}.pipes_serial PRIVATE ${PROJECT_NAME}.pipes-dispatch)

if( SIESTA_INSTALL )
  install(
    TARGETS
      ${PROJECT_NAME}.pipes_serial
      ${PROJECT_NAME}.sockets_serial
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

