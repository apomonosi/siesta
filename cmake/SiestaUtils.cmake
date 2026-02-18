
# This is borrowed from github.com/LecrisUT/CMake-Template
# Back-porting to PROJECT_IS_TOP_LEVEL to older cmake

if (CMAKE_VERSION VERSION_LESS "3.21")
  if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(PROJECT_IS_TOP_LEVEL ON)
  else ()
    set(PROJECT_IS_TOP_LEVEL OFF)
  endif ()
endif()
if (NOT DEFINED SIESTA_IS_TOP_LEVEL)
  set(SIESTA_IS_TOP_LEVEL ${PROJECT_IS_TOP_LEVEL})
endif ()

#
# Inspired by the DFTB+ project
# Stops the code if the source and the build folders are identical.
#
function(siesta_util_ensure_out_of_source_build)

  get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build this project from its source folder. Please, create a \
separate build directory and invoke CMake from that directory.")
  endif()

endfunction()


#[==========================[
Extract a certain variable to a multiline variable of certain
length.
#]==========================]
function(siesta_get_multiline)
  cmake_parse_arguments(_sgm
    "" # OPTIONS
    "LENGTH;VARIABLE;OUTPUT" # ONEVALUEARGS
    "" # MULTILINEARGS
    ${ARGN})

  if(NOT DEFINED _sgm_LENGTH)
    message(FATAL_ERROR "siesta_get_multiline must have a LENGTH argument")
  endif()

  if(NOT DEFINED _sgm_VARIABLE)
    list(POP_FRONT _sgm_UNPARSED_ARGUMENTS _sgm_VARIABLE)
  endif()
  if(NOT DEFINED _sgm_OUTPUT)
    list(POP_FRONT _sgm_UNPARSED_ARGUMENTS _sgm_OUTPUT)
  endif()


  list(JOIN ${_sgm_VARIABLE} " " _sgm_tmp)
  string(REPEAT "." ${_sgm_LENGTH} _sgm_manydots)
  string(CONCAT _sgm_manydots_regex_pattern "(" "${_sgm_manydots}" ")")
  string(REGEX REPLACE "${_sgm_manydots_regex_pattern}" "\\1\&\n\&"
    "tmp" "${_sgm_tmp}")
  set(${_sgm_OUTPUT} "${tmp}" PARENT_SCOPE)

endfunction()



#[==========================[
This function will push/pop variables to a global variable
called siesta_*_SUFFIX.
Lowercase to highlight that they are internal
This will enable one to temporarily change executable suffixes
when subsequent usages are limited.
The * should be one(or multiple) of:
- EXECUTABLE
- SHARED_LIBRARY
- STATIC_LIBRARY
#]==========================]

set(SIESTA_EXECUTABLE_SUFFIX "" CACHE STRING "Suffixes used for executables, e.g. =_omp to siesta->siesta_omp")
set(SIESTA_SHARED_LIBRARY_SUFFIX "" CACHE STRING "Suffixes used for shared libraries, e.g. =_omp to libsiesta.so->libsiesta_omp.so")
set(SIESTA_STATIC_LIBRARY_SUFFIX "" CACHE STRING "Suffixes used for static libraries, e.g. =_omp to libsiesta.a->libsiesta_omp.a")
# Internal variables (not cached)
set(siesta_EXECUTABLE_SUFFIXES "")
set(siesta_SHARED_LIBRARY_SUFFIXES "")
set(siesta_STATIC_LIBRARY_SUFFIXES "")

function(siesta_suffix)
  set(options APPEND NEW POP PREPEND)
  set(oneValueArgs SUFFIX)
  set(multiValueArgs VARIABLES)
  cmake_parse_arguments(_esop "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Get default suffix
  if(NOT DEFINED _esop_SUFFIX)
    list(LENGTH _esop_UNPARSED_ARGUMENTS nargs)
    if(nargs GREATER 0)
      list(POP_FRONT _esop_UNPARSED_ARGUMENTS _esop_SUFFIX)
    endif()
  endif()
  if(NOT DEFINED _esop_VARIABLES)
    set(_esop_VARIABLES "EXECUTABLE" "SHARED_LIBRARY" "STATIC_LIBRARY")
  endif()
  set(allowed_var "EXECUTABLE" "SHARED_LIBRARY" "STATIC_LIBRARY")
  foreach(var IN LISTS _esop_VARIABLES)
    if(NOT ${var} IN_LIST allowed_var)
      message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} expected VARIABLES to be one of ${allowed_var}, got ${var}")
    endif()
  endforeach()

  # Only one of
  #  - APPEND
  #  - NEW
  #  - POP
  #  - PREPEND
  # may be supplied
  set(defined FALSE)
  foreach(opt1 IN LISTS options)
    foreach(opt2 IN LISTS options)
      if("${opt1}" STREQUAL "${opt2}")
        if(_esop_${opt1})
          set(defined TRUE)
        endif()
      elseif(_esop_${opt1} AND _esop_${opt2})
        message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION}: both ${opt1} and ${opt2} argument found, only 1 of them allowed.")
      endif()
    endforeach()
  endforeach()
  if(NOT ${defined})
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION}: none of ${options} are defined, please at least supply one.")
  endif()


  foreach(var IN LISTS _esop_VARIABLES)

    set(cmake_var "CMAKE_${var}_SUFFIX")
    set(siesta_var "siesta_${var}_SUFFIXES")
    message(TRACE "Current ${siesta_var}=${${siesta_var}}")
    message(TRACE "Current ${cmake_var}=${${cmake_var}}")

    set(old "${${cmake_var}}")
    if(_esop_POP)
      list(POP_BACK ${siesta_var} new)
    else()
      # one of the changing of the variable
      if(_esop_APPEND)
        set(new "${old}${_esop_SUFFIX}")
      elseif(_esop_NEW)
        set(new "${_esop_SUFFIX}")
      elseif(_esop_PREPEND)
        set(new "${_esop_SUFFIX}${old}")
      endif()
      set(${siesta_var} "${${siesta_var}};${new}")
    endif()

    # Define the variables in parent scope
    set(${siesta_var} "${${siesta_var}}" PARENT_SCOPE)
    set(${cmake_var} "${new}" PARENT_SCOPE)

  endforeach()

endfunction()

macro(siesta_suffix_install)
  message(DEBUG "siesta_suffix_install running")
  siesta_suffix(
    ${SIESTA_SUFFIX}
    PREPEND
    )
  siesta_suffix(
    ${SIESTA_EXECUTABLE_SUFFIX}
    PREPEND
    VARIABLES EXECUTABLE
    )
  siesta_suffix(
    ${SIESTA_SHARED_LIBRARY_SUFFIX}
    PREPEND
    VARIABLES SHARED_LIBRARY
    )
  siesta_suffix(
    ${SIESTA_STATIC_LIBRARY_SUFFIX}
    PREPEND
    VARIABLES STATIC_LIBRARY
    )
endmacro()
macro(siesta_suffix_uninstall)
  message(DEBUG "siesta_suffix_uninstall running")
  siesta_suffix(
    ${SIESTA_SUFFIX}
    POP
    )
  siesta_suffix(
    ${SIESTA_EXECUTABLE_SUFFIX}
    POP
    VARIABLES EXECUTABLE
    )
  siesta_suffix(
    ${SIESTA_SHARED_LIBRARY_SUFFIX}
    POP
    VARIABLES SHARED_LIBRARY
    )
  siesta_suffix(
    ${SIESTA_STATIC_LIBRARY_SUFFIX}
    POP
    VARIABLES STATIC_LIBRARY
    )
endmacro()

#[==========================[
Function to dynamically add a directory as a build-dependency
if the directory exists and has some files in it
The number of files should probably be tweaked, well.

This routine accepts a set of arguments to handle how sub-directories
are added.

  DIRECTORY : required
    which directory to dynamically add
  OPTION : required
    name of the option that is exposed to the user
  NITEMS : required
    number of items that should be in DIRECTORY before it is eligeble for
    adding as a sub-project.
  HELP : optional
    help text exposed for the OPTION, defaults to a description of the OPTION and DIRECTORY
  DEFAULT : optional
    a default value of OPTION in case it is eligeble for addition.
    Defaults to TRUE.

#]==========================]
function(siesta_add_subdirectory_option)
  set(options "")
  set(oneValueArgs DIRECTORY OPTION HELP NITEMS DEFAULT)
  set(multiValueArgs "")
  cmake_parse_arguments(_asop "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT DEFINED _asop_DIRECTORY)
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} requires DIRECTORY argument")
  endif()
  if(NOT DEFINED _asop_OPTION)
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} requires OPTION argument")
  endif()
  if(NOT DEFINED _asop_NITEMS)
    set(_asop_NFILES 1)
    message(VERBOSE "${CMAKE_CURRENT_FUNCTION} set NITEMS to 1 (more than 1 file+directory in DIRECTORY)")
  endif()
  if(NOT DEFINED _asop_DEFAULT)
    set(_asop_DEFAULT TRUE)
    message(VERBOSE "${CMAKE_CURRENT_FUNCTION} set DEFAULT to TRUE if the directory exists")
  endif()
  if(NOT DEFINED _asop_HELP)
    set(_asop_HELP
      "Include support for option ${_asop_OPTION} with a default value of ${_asop_DEFAULT} "
      "if the folder ${_asop_DIRECTORY} exists with >=${_asop_NITEMS} files+directories present.")
  endif()

  message(VERBOSE "Checking directory ${_asop_DIRECTORY} for content")

  file(GLOB _result
    LIST_DIRECTORIES TRUE
    "${_asop_DIRECTORY}/*"
    )
  list(LENGTH _result _result_len)
  message(VERBOSE "Directory ${_asop_DIRECTORY} contains ${_result_len} files and directories")
  if(_result_len GREATER _asop_NITEMS)
    option(${_asop_OPTION} "${_asop_HELP}" ${_asop_DEFAULT})
    if (NOT ${_asop_DEFAULT})
      message(STATUS "Option ${_asop_OPTION} not set due to 'false' default value")
    endif()
  else()
    set(${_asop_OPTION} FALSE CACHE BOOL "${_asop_HELP}" FORCE)
  endif()
  if ( ${${_asop_OPTION}} )
    message(STATUS "Adding support with ${_asop_OPTION}=TRUE")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    add_subdirectory("${_asop_DIRECTORY}")
    list(POP_BACK CMAKE_MESSAGE_INDENT)
  endif()
endfunction()


# Create a wrapper for adding libraries and executables
# This is necessary for the Siesta build system to generalise
# the target creation.
#
# This functionality basically does:
#   target_link_libraries(name
#       ${SIESTA_LINKER_FLAGS_PRE})
#   <return action to owner>
#   target_link_libraries(name
#       ${SIESTA_LINKER_FLAGS_POST}
#       ${SIESTA_LINKER_FLAGS})
#
# However, these wrappers are just to streamline
# this process.
# This will only work in 3.19 and later due to DEFER.
#
macro(siesta_add_linker_flags)
  # Only 1 argument in this macro
  cmake_parse_arguments(_siesta_alf
    ""
    "TARGET" # single arguments
    "" # multi-key arguments
    ${ARGN})
  
  if( NOT DEFINED _siesta_alf_TARGET )
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION}: called siesta_add_linker_flags without TARGET argument")
  endif()

  # Retrieve the target name
  message(DEBUG "${CMAKE_CURRENT_FUNCTION}: Adding pre/post linker flags for ${_siesta_alf_TARGET}")
  
  list(LENGTH _siesta_alf_UNPARSED_ARGUMENTS _siesta_alf_nargs)
  if(_siesta_alf_nargs GREATER 0)
    message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION}: called siesta_add_linker_flags with too many arguments")
  endif()
  
  # Note that we default to add the PUBLIC note below.
  # This forces the definition of the target to contain
  # PUBLIC | PRIVATE | OBJECT | INTERFACE

  target_link_libraries(
    ${_siesta_alf_TARGET}
    PUBLIC
    ${SIESTA_LINKER_FLAGS_PRE}
    )

  # Create a deferred call when the current directory is left
  # Deferred calls are delayed until the final parsing.
  # This means we have to run an EVAL on the call
  cmake_language(EVAL CODE
    "cmake_language(DEFER CALL
    target_link_libraries
    ${_siesta_alf_TARGET}
    PUBLIC
    ${SIESTA_LINKER_FLAGS_POST}
    ${SIESTA_LINKER_FLAGS})"
    )

endmacro()



# Wrapper for add_executable
# Allows for easy creation of additional targets (custom/namespace).
# It also allows to add the PRE/POST linker flags implicitly.
# There are a few arguments:
#  NO_LINKER_FLAGS [option]: bypasses the addition of PRE/POST_LINKER_FLAGS
#                            If not set, will call siesta_add_linker_flags(<target>)
#  CUSTOM_TARGET <value>: defines an additional custom-target that depends on TARGET
#  NAMESPACE_TARGET <value>: defines an additional namespace target ${PROJECT_NAME}::<value> that aliases TARGET [defaults to s/${PROJECT_NAME}.\(.*\)/\1/]
#  EXPORT_NAME <value>: specifies the EXPORT_NAME property for the target [defaults to s/${PROJECT_NAME}.\(.*\)/\1/]
#  OUTPUT_NAME <value>: specifies the OUTPUT_NAME property for the target [defaults to s/${PROJECT_NAME}.\(.*\)/\1/]
function(siesta_add_executable)
  cmake_parse_arguments(_sa
    "NO_LINKER_FLAGS;NO_NAMESPACE_TARGET;NO_OUTPUT_NAME;NO_EXPORT_NAME" # options
    "TARGET;CUSTOM_TARGET;NAMESPACE_TARGET;EXPORT_NAME;OUTPUT_NAME" # single arguments
    "" # multi-key arguments
    ${ARGN})

  if( NOT DEFINED _sa_TARGET )
    list(POP_FRONT _sa_UNPARSED_ARGUMENTS _sa_TARGET)
  endif()

  # get  short target name (context need not contain PROJECT_NAME)
  string(REPLACE "${PROJECT_NAME}." "" _sa_TARGET_SHORT ${_sa_TARGET})
  
  list(APPEND CMAKE_MESSAGE_CONTEXT "${_sa_TARGET_SHORT}")
  message(DEBUG "${CMAKE_CURRENT_FUNCTION} adding target: ${_sa_TARGET}")
  list(APPEND CMAKE_MESSAGE_INDENT " - ")

  add_executable(${_sa_TARGET} ${_sa_UNPARSED_ARGUMENTS})

  # In case the IPO is enabled, do it for the target
  if(SIESTA_WITH_IPO)
    set_property(TARGET ${_sa_TARGET}
      PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()

  if( _sa_NO_LINKER_FLAGS )
    message(DEBUG "will not add linker flags to ${_sa_TARGET}")
  else()
    siesta_add_linker_flags(TARGET ${_sa_TARGET})
  endif()

  if( (NOT NO_EXPORT_NAME) AND (NOT DEFINED _sa_EXPORT_NAME) )
    set(_sa_EXPORT_NAME ${_sa_TARGET_SHORT})
  endif()
  if( DEFINED _sa_EXPORT_NAME )
    message(DEBUG "EXPORT_NAME: ${_sa_EXPORT_NAME}")
    set_target_properties(${_sa_TARGET}
      PROPERTIES
      EXPORT_NAME ${_sa_EXPORT_NAME}
      )
  endif()
  
  if( (NOT NO_OUTPUT_NAME) AND (NOT DEFINED _sa_OUTPUT_NAME) )
    set(_sa_OUTPUT_NAME ${_sa_TARGET_SHORT})
  endif()
  if( DEFINED _sa_OUTPUT_NAME )
    message(DEBUG "OUTPUT_NAME: ${_sa_OUTPUT_NAME}")
    set_target_properties(${_sa_TARGET}
      PROPERTIES
      OUTPUT_NAME ${_sa_OUTPUT_NAME}
      )
  endif()

  if( DEFINED _sa_CUSTOM_TARGET )
    message(DEBUG "CUSTOM_TARGET: ${_sa_CUSTOM_TARGET}")
    add_custom_target(${_sa_CUSTOM_TARGET} DEPENDS ${_sa_TARGET})
  endif()

  if( (NOT NO_NAMESPACE_TARGET) AND (NOT DEFINED _sa_NAMESPACE_TARGET) )
    set(_sa_NAMESPACE_TARGET ${_sa_TARGET_SHORT})
  endif()
  if( DEFINED _sa_NAMESPACE_TARGET )
    message(DEBUG "NAMESPACE_TARGET: ${PROJECT_NAME}::${_sa_NAMESPACE_TARGET}")
    add_executable(${PROJECT_NAME}::${_sa_NAMESPACE_TARGET} ALIAS ${_sa_TARGET})
  endif()
  
  list(POP_BACK CMAKE_MESSAGE_INDENT)
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)

endfunction()


# Wrapper for add_library
# Allows for easy creation of additional targets (custom/namespace).
# It also allows to add the PRE/POST linker flags implicitly.
# There are a few arguments:
#  NO_LINKER_FLAGS [option]: bypasses the addition of PRE/POST_LINKER_FLAGS
#                            If not set, will call siesta_add_linker_flags(<target>)
#  CUSTOM_TARGET <value>: defines an additional custom-target that depends on TARGET
#  NAMESPACE_TARGET <value>: defines an additional namespace target ${PROJECT_NAME}::<value> that aliases TARGET
#
# There is a an implicit addition of PIC for libraries as TYPE=OBJECT|STATIC libraries.
# This is necessary because the OBJECT|STATIC libraries are generally internal libraries
# that will be linked in the end.
# While this is not necessary for creating executables (no PIC is required)
# it is important for libsiesta and other public libraries.
# The PIC will only be added in case BUILD_SHARED_LIBS is set to TRUE.
# Otherwise the above consideration will not be applicable.
#
# NOTE:
# There is NOT an option to decide whether subsequent target_link_libraries
# requires the PUBLIC interface.
# The cmake error messages would look like this:
#
#    The keyword signature for target_link_libraries has already been used with
#    the target "<target>".  All uses of target_link_libraries with a
#    target must be either all-keyword or all-plain.
#
# This happens when one does:
#
#    target_link_libraries(<target> lib1)
#    target_link_libraries(<target> PUBLIC lib2)
#
# In this case we need to know on beforehand how to use it.
# By default we expect the 2nd invocation for all subsequent calls.
#
# A word of caution on using OBJECT libraries.
# Generally dependencies of target libraries are propagated to its dependencies.
# Hence:
#
#   add_library(foo OBJECT ...)
#   target_link_libraries(foo PUBLIC bar)
#
#   add_executable(fooexe ...)
#   target_link_libraries(fooexe PRIVATE foo)
#
# would yield different results whether it was an OBJECT or *anything else*.
# The problem arises because *object* libraries are handled as static things
# and thus INTERFACE_LIBRARIES are not propagated.
# See here:
#   https://cmake.org/cmake/help/latest/prop_tgt/INTERFACE_LINK_LIBRARIES_DIRECT.html#direct-link-dependencies-as-usage-requirements
# This is indeed very weird.
# An OBJECT library's requirements would however transition to consumers,
# and thus definitions etc. would work. However, for the OBJECT libraries
# dependencies, one would have to do something like this:
#
#   target_link_libraries(foo PUBLIC bar $<TARGET_OBJECTS:bar>)
#
# Once for the definitions etc., and once for the actual objects
# to propagate to consumers of `foo`.
# However, this becomes problematic when multiple libraries
# consumes `foo` and are then consumed in another executable:
#
#    add_executable(foobar foo another_foo_depends_on_bar)
#
# This would result in object files multiple times on the linker
# line... Yielding too many names for resolution.
function(siesta_add_library)
  cmake_parse_arguments(_sa
    "NO_LINKER_FLAGS;NO_NAMESPACE_TARGET;NO_OUTPUT_NAME;NO_EXPORT_NAME" # options
    "TARGET;CUSTOM_TARGET;NAMESPACE_TARGET;EXPORT_NAME;OUTPUT_NAME" # single arguments
    "" # multi-key arguments
    ${ARGN})

  if( NOT DEFINED _sa_TARGET )
    list(POP_FRONT _sa_UNPARSED_ARGUMENTS _sa_TARGET)
  endif()

  # get  short target name (context need not contain PROJECT_NAME)
  string(REPLACE "${PROJECT_NAME}." "" _sa_TARGET_SHORT ${_sa_TARGET})

  list(APPEND CMAKE_MESSAGE_CONTEXT "${_sa_TARGET_SHORT}")
  message(DEBUG "${CMAKE_CURRENT_FUNCTION} adding target: ${_sa_TARGET}")
  list(APPEND CMAKE_MESSAGE_INDENT " - ")

  add_library(${_sa_TARGET} ${_sa_UNPARSED_ARGUMENTS})

  # In case the IPO is enabled, do it for the target
  if(SIESTA_WITH_IPO)
    set_property(TARGET ${_sa_TARGET}
      PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()

  # Check whether we should add position independent code.
  # This is necessary since OBJECT libraries will not add this
  # implicitly...
  get_target_property(target_type ${_sa_TARGET} TYPE)
  if( SIESTA_SHARED_LIBS )
    if( "${target_type}" STREQUAL "OBJECT_LIBRARY" OR
        "${target_type}" STREQUAL "STATIC_LIBRARY") 
      set_target_properties(${_sa_TARGET} PROPERTIES
        POSITION_INDEPENDENT_CODE ${SIESTA_SHARED_LIBS}
        )
    endif()
  endif()

  if( _sa_NO_LINKER_FLAGS )
    message(DEBUG "will not add linker flags to ${_sa_TARGET}")
  else()
    siesta_add_linker_flags(TARGET ${_sa_TARGET})
  endif()

  if( (NOT NO_EXPORT_NAME) AND (NOT DEFINED _sa_EXPORT_NAME) )
    set(_sa_EXPORT_NAME ${_sa_TARGET_SHORT})
  endif()
  if( DEFINED _sa_EXPORT_NAME )
    message(DEBUG "EXPORT_NAME: ${_sa_EXPORT_NAME}")
    set_target_properties(${_sa_TARGET}
      PROPERTIES
      EXPORT_NAME ${_sa_EXPORT_NAME}
      )
  endif()
  
  if( (NOT NO_OUTPUT_NAME) AND (NOT DEFINED _sa_OUTPUT_NAME) )
    set(_sa_OUTPUT_NAME ${_sa_TARGET_SHORT})
  endif()
  if( DEFINED _sa_OUTPUT_NAME )
    message(DEBUG "OUTPUT_NAME: ${_sa_OUTPUT_NAME}")
    set_target_properties(${_sa_TARGET}
      PROPERTIES
      OUTPUT_NAME ${_sa_OUTPUT_NAME}
      )
  endif()

  if( DEFINED _sa_CUSTOM_TARGET )
    message(DEBUG "CUSTOM_TARGET: ${_sa_CUSTOM_TARGET}")
    add_custom_target(${_sa_CUSTOM_TARGET} DEPENDS ${_sa_TARGET})
  endif()

  if( (NOT NO_NAMESPACE_TARGET) AND (NOT DEFINED _sa_NAMESPACE_TARGET) )
    set(_sa_NAMESPACE_TARGET ${_sa_TARGET_SHORT})
  endif()
  if( DEFINED _sa_NAMESPACE_TARGET )
    message(DEBUG "NAMESPACE_TARGET: ${PROJECT_NAME}::${_sa_NAMESPACE_TARGET}")
    add_library(${PROJECT_NAME}::${_sa_NAMESPACE_TARGET} ALIAS ${_sa_TARGET})
  endif()
  
  list(POP_BACK CMAKE_MESSAGE_INDENT)
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)

endfunction()



# Debugging function to extract and print properties of targets
function(siesta_debug_print_properties)
  cmake_parse_arguments(_sd
    "" # options
    "TARGET" # single arguments
    "PROPERTIES" # multi-key arguments
    ${ARGN})
  
  # Get default suffix
  if(NOT DEFINED _sd_TARGET)
    list(LENGTH _sd_UNPARSED_ARGUMENTS nargs)
    if(nargs GREATER 0)
      list(POP_FRONT _sd_UNPARSED_ARGUMENTS _sd_TARGET)
    endif()
  endif()

  list(APPEND CMAKE_MESSAGE_CONTEXT "SD")
  message(DEBUG "Exploring properties for target: ${_sd_TARGET}")
  list(APPEND CMAKE_MESSAGE_INDENT "  ")
  foreach(prop IN LISTS _sd_PROPERTIES)
    get_target_property(prop_var
      ${_sd_TARGET}
      ${prop}
      )

    message(DEBUG "${prop} = ${prop_var}")

  endforeach()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  list(POP_BACK CMAKE_MESSAGE_CONTEXT)

endfunction()
