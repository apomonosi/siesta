# Print details related to the build information
# This should comprise *ALL* components.
# And preferentially in some order to retain information locally.

# This variable controls the width of the output bars
set(_pi_width 80)
set(_pi_section_package "|")
# Will allow up to 4 nested sections
set(_pi_section_delims "+" "*" "O" "@")
# create an empty list to deal with nested sections
set(_pi_section_headers)

set(_mode "STATUS")

# Global function for siesta_printing stuff
# When using the SIESTA_FIND_PACKAGE option it tells the printer
# to inform of the flags that are looked up in Siesta_find_package
# See there for details about {NAME}_FIND_METHOD, {NAME}_URL, etc.
function(siesta_print_feature_info)
  set(options REQUIRED)
  set(oneValueArgs OPTION FOUND NAME SIESTA_FIND_PACKAGE)
  set(multiValueArgs
    VARIABLES
    DEPENDENCIES OPTIONAL_DEPENDENCIES
    MSGOFF MSGON
    HEADER FOOTER
    TARGETS)
  cmake_parse_arguments(_pi "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})


  # add an empty line
  message(${_mode} "")

  if( NOT DEFINED _pi_NAME )
    message(FATAL_ERROR "NAME argument to print_feature_info is required")
  endif()

  # Default TARGETS to NAME::NAME
  if( NOT DEFINED _pi_TARGETS )
    if(TARGET "${_pi_NAME}::${_pi_NAME}")
      set(_pi_TARGETS "${_pi_NAME}::${_pi_NAME}")
    endif()
  endif()

  set(targets_exist TRUE)
  foreach(tgt ${_pi_TARGETS})
    if(NOT TARGET ${tgt})
      set(targets_exist FALSE)
      message(DEBUG "Target ${tgt} does not exist")
    else()
      message(DEBUG "Target ${tgt} exists")
    endif()
  endforeach()


  set(_pi_used "${_pi_REQUIRED}")
  set(_pi_available "${_pi_REQUIRED}")
  if(DEFINED _pi_OPTION)
    set(_pi_used ${${_pi_OPTION}})
  endif()
  if(DEFINED _pi_FOUND)
    if( ${_pi_FOUND} OR ${targets_exist})
      # Bypass the case where the user did FOUND TRUE
      # instead of FOUND SIESTA_WITH_?
      set(_pi_available TRUE)
    else()
      set(_pi_available ${${_pi_FOUND}})
    endif()
    if(NOT DEFINED _pi_OPTION)
      set(_pi_used "${_pi_available}")
    endif()
  endif()
  # streamline-value to ON/OFF
  if(_pi_used)
    set(_pi_used "ON")
  else()
    set(_pi_used "OFF")
  endif()

  siesta_print_section_line("-" 4 " ${_pi_NAME} is ${_pi_used} " ${_pi_width})
  list(APPEND CMAKE_MESSAGE_INDENT "${_pi_section_package}")

  # Add the header print-out
  if(DEFINED _pi_HEADER)
    foreach(line IN LISTS _pi_HEADER)
      message(${_mode} "${line}")
    endforeach()
  endif()
  list(APPEND CMAKE_MESSAGE_INDENT " ")

  if( _pi_used )
    if( _pi_REQUIRED )
      message(${_mode} "Required feature (cannot be disabled)")
    elseif(NOT "${_pi_OPTION}" STREQUAL "")
      message(${_mode} "Feature is turned ON and controlled by '${_pi_OPTION}'")
    endif()
  elseif( _pi_available )
    message(${_mode} "Feature is turned OFF but can be turned ON (controlled by '${_pi_OPTION}')")
  else()
    message(${_mode} "Feature is turned OFF and requires additional information to be available")
  endif()


  if( _pi_REQUIRED AND (NOT _pi_used) )
    set(_pi_tmp "")
    if( _pi_REQUIRED )
      set(_pi_tmp "${_pi_tmp}\nThis feature is a REQUIRED feature")
    endif()
    if( DEFINED _pi_OPTION)
      set(_pi_tmp "${_pi_tmp}\nFlag for controlling feature usage (is defined): ${_pi_OPTION}=${${_pi_OPTION}}")
    endif()
    if( DEFINED _pi_FOUND)
      set(_pi_tmp "${_pi_tmp}\nBuild system flag whether the package was found: ${_pi_FOUND}=${${_pi_FOUND}}")
    endif()
    message(FATAL_ERROR
      "${_pi_tmp}\n"
      "Logic in library information could not be understood."
      "The package is REQUIRED but cannot be used?")
  endif()

  # Print out information if on
  if(DEFINED _pi_MSG)
    foreach(line IN LISTS _pi_MSG)
      message(${_mode} "${line}")
    endforeach()
  endif()

  if( ${_pi_used} )
    if(DEFINED _pi_MSGON)
      foreach(line IN LISTS _pi_MSGON)
        message(${_mode} "${line}")
      endforeach()
    endif()
  elseif(DEFINED _pi_MSGOFF)
    foreach(line IN LISTS _pi_MSGOFF)
      message(${_mode} "${line}")
    endforeach()
  endif()


  if(DEFINED _pi_DEPENDENCIES)
    message(${_mode} "The following dependencies are required to use this feature:")
    foreach(dep IN LISTS _pi_DEPENDENCIES)
      message(${_mode} " - ${dep}")
    endforeach()
  endif()
  if(DEFINED _pi_OPTIONAL_DEPENDENCIES)
    message(${_mode} "The following dependencies are optional when using this feature:")
    foreach(dep IN LISTS _pi_OPTIONAL_DEPENDENCIES)
      message(${_mode} " - ${dep}")
    endforeach()
  endif()

  # Add variables if the SIESTA_FIND_PACKAGE option has been specified
  if(DEFINED _pi_SIESTA_FIND_PACKAGE)
    set(p "${_pi_SIESTA_FIND_PACKAGE}")
    list(APPEND _pi_VARIABLES
      ${p}_FOUND_METHOD
      ${p}_GIT_REPOSITORY ${p}_GIT_TAG
      ${p}_SOURCE_DIR
      )
  endif()

  macro(create_post_header header)
    macro(post_header)
      list(POP_BACK CMAKE_MESSAGE_INDENT _tmp_not_used)
      message(${_mode} "${header}")
      list(APPEND CMAKE_MESSAGE_INDENT ${_tmp_not_used})
      macro(post_header)
      endmacro()
    endmacro()
  endmacro()

  if(DEFINED _pi_VARIABLES)

    if( ${_pi_used} )

      # We have enabled the feature, so the user has supplied enough information
      create_post_header("Variables used to enable feature:")
      list(APPEND CMAKE_MESSAGE_INDENT " - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if(NOT "${${var}}" STREQUAL "")
          post_header()
          message(${_mode} "${var}=${${var}}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)

      create_post_header("Empty or undefined variables (only useful for developers!):")
      list(APPEND CMAKE_MESSAGE_INDENT " - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if("${${var}}" STREQUAL "")
          post_header()
          message(${_mode} "${var}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)

    else( ${_pi_used} )

      # We have not enabled the feature, so the user may have some missing information
      create_post_header("Variables used but the feature is NOT enabled still:")
      list(APPEND CMAKE_MESSAGE_INDENT " - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if(NOT "${${var}}" STREQUAL "")
          post_header()
          message(${_mode} "${var}=${${var}}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)

      create_post_header("Empty or undefined variables (possibly some are needed to enable the feature):")
      list(APPEND CMAKE_MESSAGE_INDENT " - ")
      foreach(var IN LISTS _pi_VARIABLES)
        if("${${var}}" STREQUAL "")
          post_header()
          message(${_mode} "${var}")
        endif()
      endforeach()
      list(POP_BACK CMAKE_MESSAGE_INDENT)

    endif( ${_pi_used} )

  endif()

  macro(_print_property_list)
    foreach(prop IN LISTS ${ARGV})
      get_target_property(out "${_target}" ${prop})
      if(NOT "${out}" STREQUAL "out-NOTFOUND")
        message(DEBUG "${prop}=${out}")
      endif()
    endforeach()
  endmacro()

  if(DEFINED _pi_TARGETS)
    message(DEBUG "Defined targets:")

    function(_print_target _target _printed_targets)

      if(NOT TARGET "${_target}")
        return()
      endif()
      if("${_target}" IN_LIST "${_printed_targets}")
        return()
      endif()

      list(APPEND "${_printed_targets}" "${_target}")
      # we need to propagate it up
      set(${_printed_targets} "${${_printed_targets}}" PARENT_SCOPE)

      message(DEBUG " - [${_target}]")
      list(APPEND CMAKE_MESSAGE_INDENT "  * ")

      # Get the type to check which variables are available
      set(_props
        NAME
        TYPE
        VERSION
        )
      _print_property_list(_props)

      get_target_property(type ${_target} TYPE)

      # global properties
      set(_props
          ALIAS_GLOBAL
          ALIASED_TARGET
          C_EXTENSIONS
          C_STANDARD
          COMPILE_DEFINITIONS
          COMPILE_FLAGS
          CXX_EXTENSIONS
          CXX_MODULE_DIRS
          CXX_STANDARD
          Fortran_MODULE_DIRECTORY
          IMPORTED
          IMPORTED_GLOBAL
          IMPORTED_LIBNAME
          IMPORTED_LINK_INTERFACE_LIBRARIES
          INCLUDE_DIRECTORIES
          INTERFACE_COMPILE_DEFINITIONS
          INTERFACE_LINK_DEPENDS
          INTERFACE_LINK_DIRECTORIES
          INTERFACE_LINK_LIBRARIES
          INTERFACE_LINK_LIBRARIES_DIRECT
          INTERFACE_LINK_OPTIONS
          INTERFACE_POSITION_INDEPENDENT_CODE
          INTERFACE_SOURCES
          INTERPROCEDURAL_OPTIMIZATION
          INSTALL_RPATH
          INSTALL_RPATH_USE_LINK_PATH
          LINK_INTERFACE_LIBRARIES
          LINK_LIBRARIES
          POSITION_INDEPENDENT_CODE
        )
      if("${type}" STREQUAL "INTERFACE_LIBRARY")

        # interface library properties
        get_target_property(_interface_libraries ${_target} INTERFACE_LINK_LIBRARIES)
        list(APPEND _props
          )

      else()

        # *standard* library properties
        get_target_property(_interface_libraries ${_target} LINK_INTERFACE_LIBRARIES)
        list(APPEND _props
          BINARY_DIR
          EXCLUDE_FROM_ALL
          IMPORTED_LIBNAME
          IMPORTED_LINK_DEPENDENT_LIBRARIES
          IMPORTED_LINK_INTERFACE_LIBRARIES
          IMPORTED_LOCATION
          IMPORTED_OBJECTS
          LINK_DEPENDS
          LINK_OPTIONS
          OUTPUT_NAME
          )

      endif()

      list(SORT _props)
      _print_property_list(_props)

      # prepare next target (nested)
      list(POP_BACK CMAKE_MESSAGE_INDENT)

      if(NOT "${_interface_libraries}" STREQUAL "_interface_libraries-NOTFOUND")
        foreach(_t IN LISTS _interface_libraries)
          _print_target("${_t}" ${_printed_targets})
        endforeach()
      endif()

    endfunction()

    # initialize a list for targets that has been printed
    set(_ptargets "")

    foreach(_target IN LISTS _pi_TARGETS)

      # Print out information for this target
      _print_target("${_target}" _ptargets)

    endforeach()
  endif()

  if(DEFINED _pi_FOOTER)
    foreach(line IN LISTS _pi_FOOTER)
      message(${_mode} "${line}")
    endforeach()
  endif()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  siesta_print_section_line("-" 0 "" ${_pi_width})

  # add an empty line
  message(${_mode} "")

endfunction()

# Section handlers

function(siesta_print_section_line delim delim_pre_count msg delim_total_count)
  # Create a message:
  #   PRE MSG POST
  # define PRE
  string(REPEAT "${delim}" ${delim_pre_count} delim_pre)
  string(LENGTH "${msg}" msg_length)
  # Calculate POST
  math(EXPR delim_post_count "${delim_total_count} - ${msg_length} - ${delim_pre_count}")
  string(REPEAT "${delim}" ${delim_post_count} delim_post)
  message(${_mode} "${delim_pre}${msg}${delim_post}")
endfunction()


macro(siesta_print_start_section _pi_msg)
  # Get current section delimiter
  list(LENGTH _pi_section_headers _pi_sec_depth)
  list(GET _pi_section_delims ${_pi_sec_depth} _pi_section_delim)

  # Append the section msg to the headers
  list(APPEND _pi_section_headers "${_pi_msg}")

  # add an empty line
  message(${_mode} "")

  # Print new section
  siesta_print_section_line("${_pi_section_delim}" 4 " ${_pi_msg} >>> " ${_pi_width})

  list(APPEND CMAKE_MESSAGE_INDENT "${_pi_section_delim} ")
endmacro()

macro(siesta_print_end_section)
  # new line, and then header
  list(POP_BACK CMAKE_MESSAGE_INDENT)

  # create the final section line
  list(POP_BACK _pi_section_headers _pi_msg)

  # Get current section delimiter
  list(LENGTH _pi_section_headers _pi_sec_depth)
  list(GET _pi_section_delims ${_pi_sec_depth} _pi_section_delim)

  # Do some arithmetic
  string(LENGTH "${_pi_msg}" _pi_section_msg_length)
  math(EXPR _pi_section_delim_len "${_pi_width} - ${_pi_section_msg_length} - 8 - 2")
  siesta_print_section_line("${_pi_section_delim}" ${_pi_section_delim_len} " <<< ${_pi_msg} " ${_pi_width})

endmacro()

macro(siesta_print_build_info)

message(${_mode} "")
message(${_mode}
  "Printing out information related to the build about to proceed. "
)
message(${_mode}
  "Please carefully go through the following lines to assert that the "
  "options are as expected, some default fall-backs may disable features "
  "depending on how variables are passed."
)

siesta_print_start_section("Compiler information")

foreach(c IN ITEMS Fortran C CXX)
  if( CMAKE_${c}_COMPILER )
    siesta_print_feature_info(
      NAME "Compiler ${c}"
      MSGON
        "Using '${CMAKE_BUILD_TYPE}' build type"
      VARIABLES
        CMAKE_${c}_COMPILER_ID
        CMAKE_${c}_COMPILER_VERSION
        CMAKE_${c}_COMPILER
        CMAKE_${c}_FLAGS
        CMAKE_${c}_FLAGS_RELEASE
        CMAKE_${c}_FLAGS_DEBUG
        CMAKE_${c}_FLAGS_MINSIZEREL
        CMAKE_${c}_FLAGS_CHECK
        CMAKE_${c}_FLAGS_RELWITHDEBINFO

      REQUIRED
     )
 endif()
endforeach()

siesta_print_feature_info(
  NAME "Linking/building"
  MSGON
    "How libraries and executables are linked"
  VARIABLES
  BUILD_SHARED_LIBS
  CMAKE_EXE_LINKER_FLAGS_INIT
  CMAKE_EXE_LINKER_FLAGS
  CMAKE_BUILD_WITH_INSTALL_RPATH
  CMAKE_INSTALL_RPATH
  CMAKE_INSTALL_RPATH_USE_LINK_PATH
  CMAKE_LIBRARY_PATH_FLAG
  CMAKE_LINK_LIBRARY_FILE_FLAG
  CMAKE_LINK_LIBRARY_FLAG
  CMAKE_MODULE_LINKER_FLAGS
  CMAKE_POSITION_INDEPENDENT_CODE
  CMAKE_SHARED_LINKER_FLAGS
  CMAKE_STATIC_LINKER_FLAGS
  CMAKE_LINK_DEPENDS_NO_SHARED
  CMAKE_SYSTEM_NAME
  CMAKE_CROSSCOMPILING

  REQUIRED
  )

siesta_print_end_section()

siesta_print_start_section("Parallel")

siesta_print_feature_info(
  NAME OpenMP
  MSGOFF
    "Can be used for extremely large systems to reduce memory requirements"
    "at the cost of some overhead since only some parts of the code is parallelized with OpenMP"
  MSGON
    "Carefully analyze your typical runs for whether this makes sense (performance wise)"
    "It might be that the overhead of using OpenMP is very high in which case it should not be used."
    "Users are encouraged to have both a non-OpenMP AND an OpenMP executable to easily switch on a case-by-case"
  VARIABLES
    SIESTA_EXECUTABLE_SUFFIX
  OPTION SIESTA_WITH_OPENMP
  )

siesta_print_feature_info(
  NAME MPI
  OPTION SIESTA_WITH_MPI
  FOUND MPI_Fortran_FOUND
  MSGOFF "Parallel support is highly advised to allow scalable and faster calculations"
  DEPENDENCIES "ScaLAPACK"
  TARGETS mpi_siesta
)

siesta_print_end_section()


siesta_print_start_section("Linear algebra")

siesta_print_feature_info(
  NAME BLAS
  VARIABLES
    BLAS_LIBRARY_DIR
    BLAS_LIBRARY
    BLAS_LIBRARIES
    BLAS_LINKER_FLAG
    BLAS_DETECTION
    BLAS_HAS_GEMM3M
  FOUND BLAS_FOUND # will generally be TRUE since otherwise the build will crash
  HEADER "Required library for fast performance"
  "Recommended libraries are:"
  "  - mkl"
  "  - openblas"
  "  - blis"
  "The NetLib BLAS library is a reference implementation which should be avoided"
  "for performance reasons"
  )

siesta_print_feature_info(REQUIRED
  NAME LAPACK
  VARIABLES
    LAPACK_LIBRARY_DIR
    LAPACK_LIBRARY
    LAPACK_LIBRARIES
    LAPACK_LINKER_FLAG
    LAPACK_DETECTION
    LAPACK_HAS_MRRR
    LAPACK_HAS_2STAGE
  FOUND LAPACK_FOUND # will generally be TRUE since otherwise the build will crash
  HEADER "Required library for fast performance"
  "Recommended libraries are:"
  "  - mkl"
  "  - openblas (can have built-in LAPACK support)"
  "  - flame"
  "The NetLib LAPACK library is fine to use as long as the linked BLAS library is"
  "NOT the NetLib BLAS library!"
  )

siesta_print_start_section("Parallel")

siesta_print_feature_info(
  NAME ScaLAPACK
  OPTION SIESTA_WITH_MPI
  FOUND SCALAPACK_FOUND
  DEPENDENCIES "MPI"
  VARIABLES
    SCALAPACK_LIBRARY_DIR
    SCALAPACK_LIBRARY
    SCALAPACK_LIBRARIES
    SCALAPACK_LINKER_FLAG
    SCALAPACK_DETECTION
    SCALAPACK_HAS_MRRR
  MSGOFF "Parallel support is highly advised to allow scalable and faster calculations"
  TARGETS SCALAPACK::SCALAPACK
  )

siesta_print_feature_info(
  NAME ELPA
  OPTION SIESTA_WITH_ELPA
  FOUND ELPA_FOUND
  MSGOFF
    "ELPA provides a significant speedup for diagonalizations, users are "
    "generally advised to add support for this library, if able."
  MSGON
    "ELPA support is controlled in fdf via:"
    "  Diag.Algorithm ELPA-2stage|ELPA-1stage # former is preferred"
  VARIABLES
    SIESTA_HAS_ELPA_THROUGH_ELSI
    SIESTA_ELPA_HAS_GPU_SUPPORT
    SIESTA_ELPA_GPU_STRING
    SIESTA_ELPA_REAL_GPU_KERNEL
    SIESTA_ELPA_COMPLEX_GPU_KERNEL
    ELPA_LIBDIR
    ELPA_LIBRARIES
    ELPA_LINK_LIBRARIES
    ELPA_INCLUDEDIR
    ELPA_INCLUDE_DIRS
    ELPA_FCFLAGS
  DEPENDENCIES "MPI" "ScaLAPACK"
  TARGETS Elpa::elpa
  )

siesta_print_feature_info(
  NAME ELSI
  OPTION SIESTA_WITH_ELSI
  FOUND ELSI_FOUND
  MSGOFF
    "ELSI provides a common interface to various solvers,"
    "including ELPA and PEXSI. Compiled on-the-fly by default."
  MSGON
    "ELSI support is controlled in fdf via:"
    "  solution-method ELSI"
    "  elsi-solver [ELPA, PEXSI, NTPoly, OMM]"
  VARIABLES
    ELSI_ELPA_HAS_GPU_SUPPORT
    ELSI_HAS_EXTERNAL_ELPA
    ELSI_HAS_PEXSI
    ELSI_ELPA_GPU_STRING
    ELSI_ELPA_REAL_GPU_KERNEL
    ELSI_ELPA_COMPLEX_GPU_KERNEL
  DEPENDENCIES "MPI" "ScaLAPACK"
  TARGETS elsi::elsi
  )

siesta_print_end_section()
siesta_print_end_section()


siesta_print_start_section("Required dependencies")

siesta_print_feature_info(REQUIRED
  NAME libfdf
  FOUND LIBFDF_FOUND
  # add VARIABLES from SIESTA_FIND_PACKAGE
  SIESTA_FIND_PACKAGE LIBFDF
  )

siesta_print_feature_info(REQUIRED
  NAME xmlf90
  FOUND XMLF90_FOUND
  SIESTA_FIND_PACKAGE XMLF90
  )

siesta_print_feature_info(REQUIRED
  NAME LibGridXC
  OPTIONAL_DEPENDENCIES Libxc MPI
  FOUND LIBGRIDXC_FOUND
  VARIABLES
    LIBGRIDXC_USES_PROCEDURE_POINTER
  SIESTA_FIND_PACKAGE LIBGRIDXC
  TARGETS libgridxc::libgridxc
  )

siesta_print_feature_info(REQUIRED
  NAME LibPSML
  HEADER "Allows using pseudo-potentials from pseudo-dojo.org"
  DEPENDENCIES xmlf90
  FOUND LIBPSML_FOUND
  VARIABLES
    LIBPSML_USES_PROCEDURE_POINTER
    LIBPSML_LINK_LIBRARIES
    LIBPSML_INCLUDE_DIRS
    LIBPSML_INCLUDEDIR
  SIESTA_FIND_PACKAGE LIBPSML
  TARGETS libpsml::libpsml
  )

siesta_print_end_section()


siesta_print_start_section("Recommended dependencies")

siesta_print_feature_info(
  NAME Libxc
  MSGOFF
    "Users are advised to use the libxc library for full feature completeness"
    "using the libxc interaction with libgridxc."
    "The XC functionals will be restricted to those directly implemented in Siesta."
  VARIABLES
    LIBXC_Fortran_INTERFACE
    LIBXC_C_LIBDIR
    LIBXC_C_LIBRARIES
    LIBXC_C_LINK_LIBRARIES
    LIBXC_C_INCLUDEDIR
    LIBXC_C_INCLUDE_DIRS
    LIBXC_F03_LIBDIR
    LIBXC_F03_LIBRARIES
    LIBXC_F03_LINK_LIBRARIES
    LIBXC_F03_INCLUDEDIR
    LIBXC_F03_INCLUDE_DIRS
    LIBXC_F90_LIBDIR
    LIBXC_F90_LIBRARIES
    LIBXC_F90_LINK_LIBRARIES
    LIBXC_F90_INCLUDEDIR
    LIBXC_F90_INCLUDE_DIRS
  OPTION SIESTA_WITH_LIBXC
  FOUND LIBXC_Fortran_FOUND
  TARGETS Libxc::xc Libxc::xc_Fortran
  )

siesta_print_feature_info(
  NAME NetCDF
  MSGOFF
    "Users are adviced to add support for NetCDF due to many functionalities"
    "relying on NetCDF file outputs."
    "It also allows parallel file outputs which can easily be compressed subsequently"
  VARIABLES
    NetCDF_ROOT
    NetCDF_PATH
    NetCDF_INCLUDE_DIR
    NetCDF_INCLUDE_DIRS
    NetCDF_LIBRARIES
    NetCDF_PARALLEL
    NetCDF_Fortran_LIBRARIES
    NetCDF_Fortran_INCLUDE_DIR
    NetCDF_Fortran_INCLUDE_DIRS

  OPTION SIESTA_WITH_NETCDF
  FOUND NetCDF_FOUND
  )

siesta_print_end_section()


siesta_print_start_section("Optional features")

siesta_print_feature_info(
  NAME flook
  MSGON "Interaction with Lua can be enabled by adding"
  "  Lua.Script luafile.lua"
  "For molecular dynamics controlled in Lua, additionally define:"
  "  MD.TypeOfRun Lua"
  VARIABLES
    FLOOK_ROOT
    FLOOK_LINK_LIBRARIES
    FLOOK_INCLUDE_DIRS
    FLOOK_SOURCE_DIR # if a submodule or fetch
  OPTION SIESTA_WITH_FLOOK
  FOUND FLOOK_FOUND
  )

siesta_print_feature_info(
  NAME pexsi
  MSGON "The native PEXSI interface can be requested by"
  "  Solution.Method PEXSI"
  VARIABLES
    PEXSI_ROOT
    PEXSI_LINK_LIBRARIES
    PEXSI_INCLUDE_DIRS
  OPTION SIESTA_WITH_PEXSI
  FOUND PEXSI_FOUND
  TARGETS ${PROJECT_NAME}::pexsi
  )

siesta_print_feature_info(
  NAME wannier90
  MSGON "The on-the-fly wrapper interface to wannier90"
  OPTION SIESTA_WITH_WANNIER90
  TARGETS wrapper-wannier90::wrapper-wannier90
  )

siesta_print_feature_info(
  NAME DFTD3
  OPTION SIESTA_WITH_DFTD3
  MSGOFF
    "Adding support for DFT-D3 can improve force descriptions"
    "using simple but heuristically determined corrections to forces"
  SIESTA_FIND_PACKAGE s-dftd3
  TARGETS s-dftd3::s-dftd3
  )


siesta_print_feature_info(
  NAME FFTW
  OPTION SIESTA_WITH_FFTW
  FOUND FFTW_DOUBLE_LIB_FOUND
  TARGETS FFTW::Double
  )


siesta_print_end_section()

siesta_print_start_section("Siesta specific options")

siesta_print_feature_info(
  NAME Suffixes
  MSGON
    "Handling the suffixes of the executables and libraries"
  VARIABLES
    SIESTA_SUFFIX
    SIESTA_SET_RPATH
    SIESTA_SHARED_LIBS
    SIESTA_EXECUTABLE_SUFFIX
    SIESTA_SHARED_LIBRARY_SUFFIX
    SIESTA_STATIC_LIBRARY_SUFFIX
  REQUIRED
  )

siesta_print_feature_info(
  NAME Kinds
  MSGON
    "Cross-compiling"
  VARIABLES
    SIESTA_WITH_HOST_OPTIMIZATION
    SIESTA_INTEGER_KINDS
    SIESTA_REAL_KINDS
  REQUIRED
  )

siesta_print_end_section()

# Empty line
message(${_mode} "Done with build information")
message(${_mode} "")

endmacro(siesta_print_build_info)
