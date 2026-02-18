# This is the main engine for SIESTA tests.
# We should probably leave only the "simple" tests as the default
# without their output verification. This can be done by choosing
# the appropriate tests when making the cmake targets (maybe via labels?)

# TODO
# Instead of doing everything in CMake, first running, then checking.
# It would be better to control the 2nd run in a wrapper executable.
# This executable would run siesta, then check.
# That would also remove the problem of labels for tests.

set(SIESTA_TESTS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/Tests")
set(SIESTA_TESTS_BINARY_DIR "${PROJECT_BINARY_DIR}/Tests")
set(SIESTA_TESTS_PSEUDO_DIR "${SIESTA_TESTS_SOURCE_DIR}/Pseudos")
set(SIESTA_UTILS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/Util")
set(SIESTA_UTILS_BINARY_DIR "${PROJECT_BINARY_DIR}/Util")


function(siesta_verifyout)
  # Parse output files.
  set(oneValueArgs NAME OUTDIR OUTFILE PARENT ISSERIAL FIXTURE)
  set(multiValueArgs LABELS)
  cmake_parse_arguments(_spout "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(DEFINED _spout_UNPARSED_ARGUMENTS)
    if(NOT DEFINED _spout_NAME)
      # get name from the first argument
      list(POP_BACK _spout_UNPARSED_ARGUMENTS _spout_NAME)
    endif()
    list(LENGTH _spout_UNPARSED_ARGUMENTS _spout_len)
    if(_spout_len GREATER 0)
      message(FATAL_ERROR "Unparsed arguments in ${CMAKE_CURRENT_FUNCTION} test=${_spout_NAME}"
        "Arguments are:\n" "  ${_spout_UNPARSED_ARGUMENTS}")
    endif()
  endif()

  set(_spout_REFDIR "${CMAKE_CURRENT_SOURCE_DIR}/Reference")
  set(name "${_spout_PARENT}[verify]")

  add_test(
    NAME ${name}
    COMMAND ${Python3_EXECUTABLE}
      "${SIESTA_TESTS_BINARY_DIR}/yaml_compare.py"
      -c "${SIESTA_TESTS_BINARY_DIR}/test-cfg.yml"
      -p "${SIESTA_TESTS_BINARY_DIR}/out_digest_yaml.awk"
      -t "${_spout_OUTDIR}/${_spout_OUTFILE}.out"
      -r "${_spout_REFDIR}/${_spout_NAME}.out"
      -o "${_spout_OUTDIR}/${_spout_NAME}.yml"
  )
  set_tests_properties(${name}
    PROPERTIES
      ENVIRONMENT "SCRIPTS_DIR='${SIESTA_TESTS_BINARY_DIR}'"
      LABELS "${_spout_LABELS};verify"
      FIXTURES_SETUP "${name}"
      FIXTURES_REQUIRED "${_spout_PARENT}"
      SKIP_RETURN_CODE 127
  )

  if(DEFINED _spout_FIXTURE)
    set(${_spout_FIXTURE} "${name}" PARENT_SCOPE)
  endif()

endfunction()


#  Function to set up tests
#
# It accepts a few arguments
# NAME : directory of test (optional, an unnamed argument will be accepted as NAME)
# PSEUDO_DIR : files to be copied in as pseudos
# EXEC: the target from which the executable would be retrieved
# MPI_NPROC : number of MPI processors
# OMP_NPROC : number of OpenMP processors
# FILES_SETUP : extra files to be copied in to the run directory
#               These files may reference files from other test as they are
#               copied in a setup step (using FIXTURES).
#               The copying step will be a FIXTURE_SETUP.
# FILES_CLEANUP : files to be deleted after running the test from the run directory
#                 The deletion step will happen as a FIXTURE_CLEANUP
# SERIAL : to skip MPI (regardless of other arguments)
# FIXTURES_SETUP|REQUIRED|CLEANUP : fixture names to assign, or require
# TEST_* : all these are return values from this script
#          It will enable one to act on the results of the script it self.
# TEST_NAME (return): name of test
# TEST_WORKING_DIRECTORY (return): working directory of test
# This function uses these external variables:
#  SIESTA_WITH_MPI
#  SIESTA_WITH_OPENMP
#  SIESTA_TESTS_MPI_NUMPROCS
### SUBTEST
function(siesta_subtest)
  set(options SERIAL)
  set(oneValueArgs NAME MPI_NPROC MPI_MINPROCS PSEUDO_DIR NAME_ALTER EXEC
    TEST_NAME
    TEST_WORKING_DIRECTORY)
  set(multiValueArgs
    FILES_SETUP FILES_CLEANUP
    LABELS EXTRA_ARGS
    FIXTURES_SETUP FIXTURES_REQUIRED FIXTURES_CLEANUP)
  cmake_parse_arguments(_stest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Check for bad arguments
  # A single unparsed-argument will be the equivalent of NAME
  # in case it has not been passed.
  if(DEFINED _stest_UNPARSED_ARGUMENTS)
    if(NOT DEFINED _stest_NAME)
      # get name from the first argument
      list(POP_BACK _stest_UNPARSED_ARGUMENTS _stest_NAME)
    endif()
    list(LENGTH _stest_UNPARSED_ARGUMENTS _stest_len)
    if(_stest_len GREATER 0)
      message(FATAL_ERROR "Unparsed arguments in ${CMAKE_CURRENT_FUNCTION} test=${_stest_NAME}"
        "Arguments are:\n"
        "  ${_stest_UNPARSED_ARGUMENTS}")
    endif()
  endif()

  # Keyword only options
  # In this case serial may be used to overwrite MPI and OpenMP tests
  if(NOT DEFINED _stest_SERIAL)
    set(_stest_SERIAL FALSE)
  endif()

  # Name of test *must* be defined
  if(NOT DEFINED _stest_NAME)
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} missing test directory argument: NAME")
  endif()

  set(msg_info "${_stest_NAME}")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  # Check for number of ranks. It will be reduced to {SIESTA_TESTS_MPI_NUMPROCS} if too high.
  if(DEFINED _stest_MPI_NPROC)
    if(SIESTA_WITH_MPI)
      if(_stest_MPI_NPROC GREATER SIESTA_TESTS_MPI_NUMPROCS)
        message(WARNING "Reducing number of MPI processors in test ${_stest_NAME} to ${SIESTA_TESTS_NUMPROCS}<${_stest_MPI_NPROC}")
        set(_stest_MPI_NPROC ${SIESTA_TESTS_MPI_NUMPROCS})
      endif()
    endif()
    set(msg_info "${msg_info} (MPI=${_stest_MPI_NPROC})")
  else()
    # set default value
    set(_stest_MPI_NPROC ${SIESTA_TESTS_MPI_NUMPROCS})
  endif()

  if(_stest_MPI_NPROC LESS _stest_MPI_MINPROCS)
    message(WARNING "Skipping test ${_stest_NAME} due to minimum number of MPI processes required.")
    return()
  endif()

  # Check for OpenMP
  # It will understand this variable as the maximum value
  if(DEFINED ENV{OMP_NUM_THREADS})
    # set the OMP_NPROC to this value if undefined
    if(DEFINED _stest_OMP_NPROC)
      if(_stest_OMP_NPROC GREATER ENV{OMP_NUM_THREADS})
        message(WARNING "Reducing number of OpenMP ranks in test ${_stest_NAME} to $ENV{OMP_NUM_THREADS}<${_stest_OMP_NPROC}")
        set(_stest_OMP_NPROC $ENV{OMP_NUM_THREADS})
      endif()
    else()
      set(_stest_OMP_NPROC $ENV{OMP_NUM_THREADS})
    endif()
    set(msg_info "${msg_info} (OMP=${_stest_OMP_NPROC})")
  else()
    # set default value of 1
    set(_stest_OMP_NPROC 1)
  endif()

  message(VERBOSE "New test ${msg_info}")

  if(NOT DEFINED _stest_PSEUDO_DIR)
    set(_stest_PSEUDO_DIR "${SIESTA_TESTS_PSEUDO_DIR}")
  else()
    message(DEBUG "test(${_stest_NAME}) will be using pseudos from: ${_stest_PSEUDO_DIR}")
  endif()

  # Build up the variables
  # First set the test name
  if(NOT DEFINED _stest_NAME_ALTER)
    set(_stest_NAME_ALTER "${_stest_NAME}")
  endif()
  set(_stest_test "${_stest_NAME_ALTER}")
  # This just does not work, using set_tests_properties does not expose the
  # env variable. Whether this is due to the command being wrapped in `sh`
  # isn't fully clear to me. I basically have to force the env in the
  # execution statement. :(
  set(_stest_env "SIESTA_PS_PATH='${_stest_PSEUDO_DIR}'") # env-vars to add
  set(_stest_pre "")
  set(_stest_post "")
  if( SIESTA_WITH_MPI AND (NOT _stest_SERIAL) AND _stest_MPI_NPROC GREATER 1 )
    set(_stest_pre
      "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_stest_MPI_NPROC} ${MPIEXEC_PREFLAGS}"
      )
    set(_stest_post
      "${MPIEXEC_POSTFLAGS}"
      )
    set(_stest_test "${_stest_test}_mpi${_stest_MPI_NPROC}")
  endif()

  if( SIESTA_WITH_OPENMP AND (NOT _stest_SERIAL) )
    list(APPEND _stest_env "OMP_NUM_THREADS=${_stest_OMP_NPROC}")
    set(_stest_test "${_stest_test}_omp${_stest_OMP_NPROC}")
  endif()

  # Create the working directory
  # Since we might have multiple tests for the same directory,
  # we will need to parse the options before doing anything
  set(_stest_wdir "${CMAKE_CURRENT_BINARY_DIR}/${_stest_test}")
  get_filename_component(_stest_parent_dir "${CMAKE_CURRENT_BINARY_DIR}" NAME)

  # Allow for user defined executables (through target names)
  if(NOT DEFINED _stest_EXEC)
    set(_stest_EXEC "${PROJECT_NAME}.siesta")
  endif()
  get_target_property(_stest_EXEC_name ${_stest_EXEC} OUTPUT_NAME)
  if( "${_stest_EXEC_name}" STREQUAL "" )
    set(_stest_EXEC_name "${_stest_EXEC}")
  endif()
  set(_stest_FULL_NAME "${_stest_EXEC_name}-${_stest_parent_dir}-${_stest_test}")


  # Define the pre-creation steps

  # This test creates the directories, and adds it to the FULL_NAME-setup fixture
  # It creates the directories for the test to run. So only at runtime will it be
  # created.
  # In this case, the directories are full paths, so there is no requirement for
  # a working directory clause!
  add_test(
    NAME "${_stest_FULL_NAME}-mkdir"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_stest_wdir}"
    COMMAND_EXPAND_LISTS
  )
  set_tests_properties("${_stest_FULL_NAME}-mkdir"
    PROPERTIES
      # Add it to the -setup fixture
      FIXTURES_SETUP "${_stest_FULL_NAME}-setup"
  )

  # Now we can prepare the test directory
  # Get the list of fdf-files
  file(GLOB fdf_files "*.fdf")
  message(DEBUG "test(${_stest_NAME}) copying fdf files: ${fdf_files}")
  # Check if the main fdf file is found.
  # If not, this test cannot be created.
  list(FIND fdf_files "${CMAKE_CURRENT_SOURCE_DIR}/${_stest_NAME}.fdf" main_idx)
  if( ${main_idx} LESS 0 )
    message(FATAL_ERROR
      "Could not create test for ${_stest_NAME} since the main "
      "fdf file could not be found."
      ""
      "Searched for ${_stest_NAME}.fdf in '${CMAKE_CURRENT_SOURCE_DIR}'"
    )
  endif()

  # Since we are copying files from a specific source-directory,
  # we have to explicitly define the source-dir
  add_test(
    NAME "${_stest_FULL_NAME}-fdf-setup"
    COMMAND ${CMAKE_COMMAND} -E copy ${fdf_files} "${CMAKE_CURRENT_BINARY_DIR}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMAND_EXPAND_LISTS # expands list of files if they are lists
  )
  set_tests_properties("${_stest_FULL_NAME}-fdf-setup"
    PROPERTIES
      # We can't copy before the directory has been created.
      # So it depends on that one, but they are both part of the
      # same -setup fixture
      # Multiple tests can belong to the same fixture name, in which
      # case they can run in any order, the DEPEND clause specifies order.
      DEPENDS "${_stest_FULL_NAME}-mkdir"
      # Add it to the -setup fixture
      FIXTURES_SETUP "${_stest_FULL_NAME}-setup"
  )

  if(DEFINED _stest_FILES_SETUP)
    message(DEBUG "test(${_stest_NAME}) copying setup files: ${_stest_FILES_SETUP}")
    add_test(
      NAME "${_stest_FULL_NAME}-files-setup"
      COMMAND ${CMAKE_COMMAND} -E copy ${_stest_FILES_SETUP} "${_stest_wdir}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      COMMAND_EXPAND_LISTS # expands list of files if they are lists
    )
    set_tests_properties("${_stest_FULL_NAME}-files-setup"
      PROPERTIES
        DEPENDS "${_stest_FULL_NAME}-mkdir"
        # Add it to the -setup fixture
        FIXTURES_SETUP "${_stest_FULL_NAME}-setup"
    )
  endif()

  # Cleanup files after run
  if(DEFINED _stest_FILES_CLEANUP)
    message(DEBUG "test(${_stest_NAME}) removing files: ${_stest_FILES_CLEANUP}")
    add_test(
      NAME "${_stest_FULL_NAME}-files-cleanup"
      COMMAND ${CMAKE_COMMAND} -E rm -f ${_stest_FILES_CLEANUP}
      WORKING_DIRECTORY "${_stest_wdir}"
      COMMAND_EXPAND_LISTS # expands files if more are present
    )
    set_tests_properties("${_stest_FULL_NAME}-files-cleanup"
      PROPERTIES
        # Add it to the -cleanup fixture
        FIXTURES_CLEANUP "${_stest_FULL_NAME}-cleanup"
    )
  endif()

  add_test(
    NAME "${_stest_FULL_NAME}"
    COMMAND sh -c "SIESTA_PS_PATH='${_stest_PSEUDO_DIR}' ${_stest_pre} \
        '$<TARGET_FILE:${_stest_EXEC}>' -out ${_stest_test}.out ${_stest_EXTRA_ARGS} \
        ../${_stest_NAME}.fdf ${_stest_post}"
    WORKING_DIRECTORY "${_stest_wdir}"
    COMMAND_EXPAND_LISTS
  )

  # Add environment variables to the test
  set(fixtures)
  # Adding them manually in one big list failed to include
  # them correctly, so now we create the list of fixtures,
  # so we add them correctly.
  list(APPEND fixtures ${_stest_FIXTURES_SETUP})
  list(APPEND fixtures ${_stest_FIXTURES_CLEANUP})
  list(APPEND fixtures ${_stest_FIXTURES_REQUIRES})
  # This tests setup action, mkdir -> copy fdf & copy FILES
  list(APPEND fixtures ${_stest_FULL_NAME}-setup)
  # This tests cleanup action, rm files
  list(APPEND fixtures ${_stest_FULL_NAME}-cleanup)
  set_tests_properties("${_stest_FULL_NAME}"
    PROPERTIES
      ENVIRONMENT "${_stest_env}"
      LABELS "${_stest_LABELS}"
      # We create this test as a SETUP fixture as well
      # This will enable the parse-out to use this as a fixture.
      # Recall, that depend does not check whether the test succeeded or not.
      # A DEPEND will after the dependency has executed, regardless of success.
      FIXTURES_SETUP "${_stest_FULL_NAME}"
      FIXTURES_REQUIRED "${fixtures}"
      SKIP_RETURN_CODE 127
  )

  # Add output verification
  siesta_verifyout(${_stest_NAME_ALTER}
    ISSERIAL ${_stest_SERIAL}
    # Note, the PARENT name is also the name of the fixture
    # that is being REQUIRED in siesta_verifyout
    PARENT "${_stest_FULL_NAME}"
    OUTFILE "${_stest_test}"
    OUTDIR "${_stest_wdir}"
    LABELS "${_stest_LABELS}"
  )

  if(DEFINED _stest_TEST_NAME)
    set(${_stest_TEST_NAME} "${_stest_FULL_NAME}" PARENT_SCOPE)
  endif()
  if(DEFINED _stest_TEST_WORKING_DIRECTORY)
    set(${_stest_TEST_WORKING_DIRECTORY} "${_stest_wdir}" PARENT_SCOPE)
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)

endfunction()

# Function to test utilities
function(siesta_util_test)
  set(options)
  set(oneValueArgs NAME BASH_SCR BASH_VERIF EXECUTABLE INPUT_FILE EXEC_DIR
                   SOURCE_DIR OUTFILE REFDIR)
  set(multiValueArgs LABELS)
  cmake_parse_arguments(_utils "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Check for bad arguments
  if(DEFINED _utils_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unparsed arguments in ${CMAKE_CURRENT_FUNCTION} test=${_utils_NAME}"
      "Arguments are:\n"
      "  ${_utils_UNPARSED_ARGUMENTS}")
  endif()

  # Name of test *must* be defined
  if(NOT DEFINED _utils_NAME)
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} missing test name argument")
  endif()

  # Set the input directory
  set(_utils_inputs_dir "${_utils_SOURCE_DIR}")

  # Retrieve the name of the input file and bash script
  get_filename_component(_utils_INPUT_FILENAME ${_utils_INPUT_FILE} NAME)
  get_filename_component(_utils_bash_test ${_utils_BASH_SCR} NAME)

  # Create a unique test name
  set(_utils_FULL_NAME "Utility_Test_${_utils_NAME}")

  # Create the working directory
  set(_utils_wdir "${_utils_EXEC_DIR}/Tests/${_utils_FULL_NAME}")

  # Create mkdir test
  add_test(
    NAME "${_utils_FULL_NAME}-mkdir"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_utils_wdir}"
  )
  set_tests_properties("${_utils_FULL_NAME}-mkdir"
    PROPERTIES
      FIXTURES_SETUP "${_utils_FULL_NAME}-setup"
  )

  # Copy the input file and the bash script to the build directory for the specific utility test
  configure_file(${_utils_INPUT_FILE} ${_utils_wdir}/${_utils_INPUT_FILENAME} COPYONLY)
  configure_file(${_utils_BASH_SCR} ${_utils_wdir}/${_utils_bash_test} COPYONLY)

  # Create the actual utility test
  add_test(
    NAME "${_utils_FULL_NAME}"
    COMMAND sh ${_utils_wdir}/${_utils_bash_test}
      ${_utils_EXECUTABLE}
      ${_utils_INPUT_FILENAME}
      ${_utils_wdir}
      ${_utils_SOURCE_DIR}
      ${_utils_EXEC_DIR}
      ${_utils_OUTFILE}
    WORKING_DIRECTORY "${_utils_wdir}"
  )

  set_tests_properties("${_utils_FULL_NAME}"
    PROPERTIES
      FIXTURES_SETUP "${_utils_FULL_NAME}"
      FIXTURES_REQUIRED "${_utils_FULL_NAME}-setup"
      LABELS "${_utils_LABELS}"
  )

  if (DEFINED _utils_BASH_VERIF)
    get_filename_component(_utils_bash_verifier ${_utils_BASH_VERIF} NAME)
    configure_file(${_utils_BASH_VERIF} ${_utils_wdir}/${_utils_bash_verifier} COPYONLY)

    add_test(
      NAME "${_utils_FULL_NAME}_verify"
      COMMAND sh ${_utils_wdir}/${_utils_bash_verifier}
      ${_utils_EXECUTABLE}
      ${_utils_wdir}
      ${_utils_OUTFILE}
      ${_utils_REFDIR}
      WORKING_DIRECTORY "${_utils_wdir}"
    )

    set_tests_properties("${_utils_FULL_NAME}_verify"
    PROPERTIES
      FIXTURES_SETUP "${_utils_FULL_NAME}-verify"
      FIXTURES_REQUIRED "${_utils_FULL_NAME}"
      LABELS "${_utils_LABELS}"
      SKIP_RETURN_CODE 127
    )
  endif()
endfunction()
#  Function to set up tests, similar to siesta_subtest.
#
# It accepts a few arguments
# NAME : directory of test (optional, an unnamed argument will be accepted as NAME)
# PSEUDO_DIR : files to be copied in as pseudos
# EXEC: the target from which the executable would be retrieved
# EXEC_DRIVER: the target from which the executable qmmm driver would be retrieved
# MPI_NPROC : number of MPI processors
# OMP_NPROC : number of OpenMP processors
# FILES_SETUP : extra files to be copied in to the run directory
# SERIAL : to skip MPI (regardless of other arguments)
function(siestaqmmm_subtest)
  set(options
      SERIAL MMONLY)
  set(oneValueArgs
      NAME MPI_NPROC PSEUDO_DIR NAME_ALTER
      EXEC EXEC_DRIVER
      TEST_NAME TEST_WORKING_DIRECTORY)
  set(multiValueArgs
      FILES_SETUP FILES_CLEANUP
      LABELS EXTRA_ARGS
      FIXTURES_SETUP FIXTURES_REQUIRED FIXTURES_CLEANUP)
  cmake_parse_arguments(_stest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Check for bad arguments
  # A single unparsed-argument will be the equivalent of NAME
  # in case it has not been passed.
  if(DEFINED _stest_UNPARSED_ARGUMENTS)
    if(NOT DEFINED _stest_NAME)
      # get name from the first argument
      list(POP_BACK _stest_UNPARSED_ARGUMENTS _stest_NAME)
    endif()
    list(LENGTH _stest_UNPARSED_ARGUMENTS _stest_len)
    if(_stest_len GREATER 0)
      message(FATAL_ERROR "Unparsed arguments in ${CMAKE_CURRENT_FUNCTION} test=${_stest_NAME}"
        "Arguments are:\n"
        "  ${_stest_UNPARSED_ARGUMENTS}")
    endif()
  endif()

  # Keyword only options
  # In this case serial may be used to overwrite MPI and OpenMP tests
  if(NOT DEFINED _stest_SERIAL)
    set(_stest_SERIAL FALSE)
  endif()

  if(NOT DEFINED _stest_MMONLY)
    set(_stest_MMONLY FALSE)
  endif()

  # Name of test *must* be defined
  if(NOT DEFINED _stest_NAME)
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} missing test directory argument: NAME")
  endif()

  set(msg_info "${_stest_NAME}")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")

  # Check for number of ranks. It will be reduced to {SIESTA_TESTS_MPI_NUMPROCS} if too high.
  if(DEFINED _stest_MPI_NPROC)
    if(SIESTA_WITH_MPI)
      if(_stest_MPI_NPROC GREATER SIESTA_TESTS_MPI_NUMPROCS)
        message(WARNING "Reducing number of MPI processors in test ${_stest_NAME} to ${SIESTA_TESTS_NUMPROCS}<${_stest_MPI_NPROC}")
        set(_stest_MPI_NPROC ${SIESTA_TESTS_MPI_NUMPROCS})
      endif()
    endif()
    set(msg_info "${msg_info} (MPI=${_stest_MPI_NPROC})")
  else()
    # set default value
    set(_stest_MPI_NPROC ${SIESTA_TESTS_MPI_NUMPROCS})
  endif()

  # Check for OpenMP
  # It will understand this variable as the maximum value
  if(DEFINED ENV{OMP_NUM_THREADS})
    # set the OMP_NPROC to this value if undefined
    if(DEFINED _stest_OMP_NPROC)
      if(_stest_OMP_NPROC GREATER ENV{OMP_NUM_THREADS})
        message(WARNING "Reducing number of OpenMP ranks in test ${_stest_NAME} to $ENV{OMP_NUM_THREADS}<${_stest_OMP_NPROC}")
        set(_stest_OMP_NPROC $ENV{OMP_NUM_THREADS})
      endif()
    else()
      set(_stest_OMP_NPROC $ENV{OMP_NUM_THREADS})
    endif()
    set(msg_info "${msg_info} (OMP=${_stest_OMP_NPROC})")
  else()
    # set default value of 1
    set(_stest_OMP_NPROC 1)
  endif()

  message(VERBOSE "New test ${msg_info}")

  if(NOT DEFINED _stest_PSEUDO_DIR)
    set(_stest_PSEUDO_DIR "${SIESTA_TESTS_PSEUDO_DIR}")
  else()
    message(DEBUG "test(${_stest_NAME}) will be using pseudos from: ${_stest_PSEUDO_DIR}")
  endif()

  # Build up the variables
  # First set the test name
  if(NOT DEFINED _stest_NAME_ALTER)
    set(_stest_NAME_ALTER "${_stest_NAME}")
  endif()
  set(_stest_test "${_stest_NAME_ALTER}")
  # This just does not work, using set_tests_properties does not expose the
  # env variable. Whether this is due to the command being wrapped in `sh`
  # isn't fully clear to me. I basically have to force the env in the
  # execution statement. :(
  set(_stest_env "SIESTA_PS_PATH='${_stest_PSEUDO_DIR}'") # env-vars to add
  set(_stest_pre "")
  set(_stest_post "")
  if( SIESTA_WITH_MPI AND (NOT _stest_SERIAL) AND _stest_MPI_NPROC GREATER 1 )
    set(_stest_pre
      "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_stest_MPI_NPROC} ${MPIEXEC_PREFLAGS}"
      )
    set(_stest_post
      "${MPIEXEC_POSTFLAGS}"
      )
    set(_stest_test "${_stest_test}_mpi${_stest_MPI_NPROC}")
  endif()

  if( SIESTA_WITH_OPENMP AND (NOT _stest_SERIAL) )
    list(APPEND _stest_env "OMP_NUM_THREADS=${_stest_OMP_NPROC}")
    set(_stest_test "${_stest_test}_omp${_stest_OMP_NPROC}")
  endif()

  # Create the working directory
  # Since we might have multiple tests for the same directory,
  # we will need to parse the options before doing anything
  set(_stest_wdir "${CMAKE_CURRENT_BINARY_DIR}/${_stest_test}")
  get_filename_component(_stest_parent_dir "${CMAKE_CURRENT_BINARY_DIR}" NAME)

  # Dirty fix for transiesta test.
  if(NOT DEFINED _stest_EXEC)
    set(_stest_EXEC "${PROJECT_NAME}.siesta")
  endif()
  if(NOT DEFINED _stest_EXEC_DRIVER)
    set(_stest_EXEC_DRIVER "${PROJECT_NAME}.siesta_qmmm")
  endif()

  get_target_property(_stest_EXEC_name ${_stest_EXEC_DRIVER} OUTPUT_NAME)
  if( "${_stest_EXEC_name}" STREQUAL "" )
    set(_stest_EXEC_name "${_stest_EXEC_DRIVER}")
  endif()
  set(_stest_FULL_NAME "${_stest_EXEC_name}-${_stest_parent_dir}-${_stest_test}")

  # Now we can prepare the test directory:
  add_test(
    NAME "${_stest_FULL_NAME}-mkdir"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_stest_wdir}"
    COMMAND_EXPAND_LISTS
  )
  set_tests_properties("${_stest_FULL_NAME}-mkdir"
    PROPERTIES
      # Add it to the -setup fixture
      FIXTURES_SETUP "${_stest_FULL_NAME}-setup"
  )

  file(GLOB fdf_files "*.fdf")
  file(COPY
    ${fdf_files}
    ../testscript-qmmm.sh
    ../testscript-mm.sh
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
    )

  # Since we are copying files from a specific source-directory,
  # we have to explicitly define the source-dir
  add_test(
    NAME "${_stest_FULL_NAME}-fdf-setup"
    COMMAND ${CMAKE_COMMAND} -E copy ${fdf_files} "${CMAKE_CURRENT_BINARY_DIR}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMAND_EXPAND_LISTS # expands list of files if they are lists
  )
  set_tests_properties("${_stest_FULL_NAME}-fdf-setup"
    PROPERTIES
      # We can't copy before the directory has been created.
      # So it depends on that one, but they are both part of the
      # same -setup fixture
      # Multiple tests can belong to the same fixture name, in which
      # case they can run in any order, the DEPEND clause specifies order.
      DEPENDS "${_stest_FULL_NAME}-mkdir"
      # Add it to the -setup fixture
      FIXTURES_SETUP "${_stest_FULL_NAME}-setup"
  )

  if(DEFINED _stest_FILES_SETUP)
  message(DEBUG "test(${_stest_FULL_NAME}) copying setup files: ${_stest_FILES_SETUP}")
  add_test(
    NAME "${_stest_FULL_NAME}-files-setup"
    COMMAND ${CMAKE_COMMAND} -E copy ${_stest_FILES_SETUP} "${_stest_wdir}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMAND_EXPAND_LISTS # expands list of files if they are lists
  )
  set_tests_properties("${_stest_FULL_NAME}-files-setup"
    PROPERTIES
      DEPENDS "${_stest_FULL_NAME}-mkdir"
      # Add it to the -setup fixture
      FIXTURES_SETUP "${_stest_FULL_NAME}-setup"
  )
  endif()

  # Cleanup files after run
  if(DEFINED _stest_FILES_CLEANUP)
    message(DEBUG "test(${_stest_NAME}) removing files: ${_stest_FILES_CLEANUP}")
    add_test(
      NAME "${_stest_FULL_NAME}-files-cleanup"
      COMMAND ${CMAKE_COMMAND} -E rm -f ${_stest_FILES_CLEANUP}
      WORKING_DIRECTORY "${_stest_wdir}"
      COMMAND_EXPAND_LISTS # expands files if more are present
    )
    set_tests_properties("${_stest_FULL_NAME}-files-cleanup"
      PROPERTIES
        # Add it to the -cleanup fixture
        FIXTURES_CLEANUP "${_stest_FULL_NAME}-cleanup"
    )
  endif()

  if(NOT DEFINED _stest_SCRIPT)
    if ( NOT _stest_MMONLY )
      set(_stest_SCRIPT "../testscript-qmmm.sh")
    else()
      set(_stest_SCRIPT "../testscript-mm.sh")
    endif()
  endif()

  if ( NOT _stest_MMONLY )
    add_test(
      NAME "${_stest_FULL_NAME}"
      COMMAND sh -c "${_stest_SCRIPT} '$<TARGET_FILE:${_stest_EXEC_DRIVER}>' \
        '../${_stest_NAME}.fdf' '${_stest_test}.out' '${_stest_PSEUDO_DIR}' \
        '$<TARGET_FILE:${_stest_EXEC}>' '../${_stest_NAME}.siestain.fdf' \
        '${_stest_test}.siesta.out' '${_stest_pre}' '${_stest_post}' "
      WORKING_DIRECTORY "${_stest_wdir}"  )
  else()
    add_test(
      NAME "${_stest_FULL_NAME}"
      COMMAND sh -c "${_stest_SCRIPT} '$<TARGET_FILE:${_stest_EXEC_DRIVER}>' \
        '../${_stest_NAME}.fdf' '${_stest_test}.out' "
      WORKING_DIRECTORY "${_stest_wdir}"  )
  endif()

  # Add environment variables to the test
  set(fixtures)
  list(APPEND fixtures ${_stest_FIXTURES_SETUP})
  list(APPEND fixtures ${_stest_FIXTURES_CLEANUP})
  list(APPEND fixtures ${_stest_FIXTURES_REQUIRES})
  list(APPEND fixtures ${_stest_FULL_NAME}-setup)
  list(APPEND fixtures ${_stest_FULL_NAME}-cleanup)

  set_tests_properties("${_stest_FULL_NAME}"
    PROPERTIES
      ENVIRONMENT "${_stest_env}"
      LABELS "${_stest_LABELS}"
      FIXTURES_SETUP "${_stest_FULL_NAME}"
      FIXTURES_REQUIRED "${fixtures}"
      SKIP_RETURN_CODE 127)

  # Possible way to enforce dependency: FIXTURES_SETUP "${_stest_test}-dep"
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  # Add output verification
  siesta_verifyout(${_stest_NAME_ALTER}
    ISSERIAL ${_stest_SERIAL}
    PARENT "${_stest_FULL_NAME}"
    OUTFILE "${_stest_test}"
    OUTDIR "${_stest_wdir}"
    LABELS "${_stest_LABELS}"
    )

endfunction()
