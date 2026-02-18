
set(SIESTA_TESTS_MPI_MIN_NUMPROCS "1" CACHE STRING "Minimum number of ranks used for MPI tests")
if( SIESTA_WITH_MPI )
  set(SIESTA_TESTS_MPI_MAX_NUMPROCS "${MPIEXEC_MAX_NUMPROCS}" CACHE STRING "Maximum number of ranks used for MPI tests")
else()
  set(SIESTA_TESTS_MPI_MAX_NUMPROCS "1" CACHE STRING "Maximum number of ranks used for MPI tests")
endif()

# NOT CACHE! Will this even work when ctest is runned a 2nd time?
# I.e. won't all options be cached at config time?
# Allows one to rerun it
if( SIESTA_TESTS_MPI_MAX_NUMPROCS GREATER 4 )
  set(SIESTA_TESTS_MPI_NUMPROCS "4")
else()
  set(SIESTA_TESTS_MPI_NUMPROCS "${SIESTA_TESTS_MPI_MAX_NUMPROCS}")
endif()


# This function uses these external variables:
# SIESTA_TESTS_MPI_MIN_NUMPROCS
# SIESTA_TESTS_MPI_NUMPROCS
# SIESTA_TESTS_MPI_MAX_NUMPROCS
function(siesta_test_get_mpi)
  # Define options that are exclusive
  # These are just wrappers for getting MPI variables (if defined)
  # The CMD_LINE variable will be a command line that combines these:
  #  MPIEXEC_EXECUTABLE
  #    MPIEXEC_PREFLAGS
  #    MPIEXEC_NUMPROC_FLAG <numprocs>
  #    EXECUTABLE
  #    MPIEXEC_POSTFLAGS
  #    FLAGS
  set(options_get
    MPIEXEC_EXECUTABLE
    MPIEXEC_PREFLAGS
    MPIEXEC_NUMPROC_FLAG
    MPIEXEC_MAX_NUMPROCS
    NUMPROCS # number of processors used (corrected for the MAX_NUMPROCS and otherwise)
    COMMAND # corrected COMMAND line
    MPIEXEC_POSTFLAGS
    )
    
  set(options "${options_get}")
  set(oneValueArgs
    OUTPUT
    MIN_NUMPROCS
    MAX_NUMPROCS
    EXECUTABLE
    )
  set(multiValueArgs FLAGS)

  cmake_parse_arguments(_stmpi "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Default OUTPUT to the unnamed argument
  if( DEFINED _stmpi_UNPARSED_ARGUMENTS )
    if( NOT DEFINED _stmpi_OUTPUT )
      list(POP_BACK _stmpi_OUTPUT _stmpi_UNPARSED_ARGUMENTS)
    endif()
  
    list(LENGTH _stmpi_UNPARSED_ARGUMENTS _stmpi_argv_len)
    if( _stmpi_argv_len GREATER 0 )
      message(FATAL_ERROR "Unparsed arguments in ${CMAKE_CURRENT_FUNCTION}: "
        "Arguments are:\n" "  ${_stmpi_UNPARSED_ARGUMENTS}")
    endif()
  endif()

  if( NOT DEFINED _stmpi_OUTPUT )
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} did not define an output variable")
  endif()
  
  # Populate defaults
  if( NOT DEFINED _stmpi_MIN_NUMPROCS )
    set(_stmpi_MIN_NUMPROCS "${SIESTA_TESTS_MPI_MIN_NUMPROCS}")
  endif()
  if( NOT DEFINED _stmpi_MAX_NUMPROCS )
    set(_stmpi_MAX_NUMPROCS "${SIESTA_TESTS_MPI_MAX_NUMPROCS}")
  endif()

  # TODO check for MIN/MAX_NUMPROCS against global variables SIESTA_TESTS_MPI_*_NUMPROCS
  set(numprocs "${SIESTA_TESTS_MPI_NUMPROCS}")

  # Define the actual number of processors we can use
  if( ${SIESTA_TESTS_MPI_NUMPROCS} GREATER ${_stmpi_MAX_NUMPROCS})
    # This corner case will simply truncate the allowed value.
    # if there is no consistency between MIN/MAX, then the MIN test will fail.
    set(numprocs "${_stmpi_MAX_NUMPROCS}")
  endif()
  if( ${SIESTA_TESTS_MPI_NUMPROCS} LESS ${_stmpi_MIN_NUMPROCS})
    # TODO, this corner case should be handled differently, perhaps
    # an ELIGIBLE request, T|F to quickly check if it is valid, or not?
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} requested fewer resources than the minimum allowed.")
  endif()

  # Ensure exclusive usage
  set(i 0)
  foreach(name IN ITEMS ${options_get})
    if( _stmpi_${name} )
      math(EXPR i "${i}+1")
    endif()
  endforeach()
  if( NOT ${i} EQUAL 1 )
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} got either 0 or more than 1 MPI requests, only use 1!")
  endif()

  if( NOT SIESTA_WITH_MPI )
    # Quick return
    if( _stmpi_NUMPROCS OR _stmpi_MPIEXEC_MAX_NUMPROCS)
      set(${_stmpi_OUTPUT} 1 PARENT_SCOPE)
    elseif( _stmpi_COMMAND )
      set(${_stmpi_OUTPUT} "${_stmpi_EXECUTABLE} ${_stmpi_FLAGS}" PARENT_SCOPE)
    else()
      set(${_stmpi_OUTPUT} "" PARENT_SCOPE)
    endif()
    return()
  endif()

  # Check the ranks

  # We are in an MPI environment, and now we need to do by
  if( _stmpi_MPIEXEC_EXECUTABLE )
    set(${_stmpi_OUTPUT} "${MPIEXEC_EXECUTABLE}" PARENT_SCOPE)
  elseif( _stmpi_MPIEXEC_PREFLAGS )
    set(${_stmpi_OUTPUT} "${MPIEXEC_PREFLAGS}" PARENT_SCOPE)
  elseif( _stmpi_MPIEXEC_NUMPROC_FLAG )
    set(${_stmpi_OUTPUT} "${MPIEXEC_NUMPROC_FLAG}" PARENT_SCOPE)
  elseif( _stmpi_NUMPROCS )
    set(${_stmpi_OUTPUT} "${numprocs}" PARENT_SCOPE)
  elseif( _stmpi_MPIEXEC_MAX_NUMPROCS )
    set(${_stmpi_OUTPUT} "${MPIEXEC_MAX_NUMPROCS}" PARENT_SCOPE)
  elseif( _stmpi_COMMAND )
    if( DEFINED _stmpi_EXECUTABLE )
      set(${_stmpi_OUTPUT} "${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPI_NUMPROC_FLAG} ${numprocs} ${_stmpi_EXECUTABLE} ${MPIEXEC_POSTFLAGS} ${_stmpi_ARGS}" PARENT_SCOPE)
    else()
      message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} cannot return a command line without an EXECUTABLE!")
    endif()
  elseif( _stmpi_MPIEXEC_POSTFLAGS )
    set(${_stmpi_OUTPUT} "${MPIEXEC_POSTFLAGS}" PARENT_SCOPE)
  else()
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} should never happen!")
  endif()

endfunction()
