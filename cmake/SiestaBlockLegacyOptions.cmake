#------------------------------------------------------
#
#  This file contains interim functionality to detect the legacy form of options
#  and other variables (e.g. WITH_MPI instead of the preferred SIESTA_WITH_MPI).

message(CHECK_START "Check use of deprecated variables")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

message(STATUS "All Siesta variables/options are prefixed with SIESTA_, e.g. SIESTA_WITH_MPI")

set(_dep_var_succeed TRUE)

foreach(opt
	WITH_MPI
	WITH_OPENMP
	WITH_NETCDF
	WITH_NCDF
	WITH_NETCDF_PARALLEL
	WITH_LIBXC
	WITH_GRID_SP
	WITH_PEXSI
	WITH_DFTD3
	WITH_DFTD3_TESTS
	WITH_WANNIER90
	WITH_CHESS
	WITH_ELPA
	WITH_ELSI
	WITH_FLOOK
	WITH_FFTW
	WITH_NO_MPI_INTERFACES
	WITH_PROFILE_NVTX
	BUILD_DOCS
)

  if( DEFINED ${opt} )
    if ( NOT DEFINED SIESTA_${opt})
      if (PROJECT_IS_TOP_LEVEL)
        message(FATAL_ERROR "Using deprecated ${opt} option variable instead of SIESTA_${opt}")
      else()
        message(WARNING "${opt} is a deprecated option name in Siesta. Define SIESTA_${opt} if needed")
        set(_dep_var_succeed FALSE)
      endif()
    elseif(NOT "${${opt}}" STREQUAL "${SIESTA_${opt}}")

      message(STATUS "- Both the deprecated ${opt} option variable name and the new variable SIESTA_${opt} are set")
      message(STATUS "  SIESTA_${opt}: ${SIESTA_${opt}}")
      message(STATUS "  ${opt}: ${${opt}}")
      message(WARNING "The new variable SIESTA_${opt} takes precedence")
      set(_dep_var_succeed FALSE)

    endif()
  endif()

endforeach()

foreach(var
  WITH_UNIT_CONVENTION
	PROFILE_NVTX_LIBRARY
)

  if( DEFINED ${var} )
    if ( NOT DEFINED SIESTA_${var})

      if (PROJECT_IS_TOP_LEVEL)
        message(FATAL_ERROR "Using deprecated ${var} variable name instead of SIESTA_${var}")
      else()
        message(WARNING "${var} is a deprecated variable name in Siesta. Define SIESTA_${var} if needed")
        set(_dep_var_succeed FALSE)
      endif()

    elseif(NOT "${${var}}" STREQUAL "${SIESTA_${var}}")

      message(STATUS " >> Both the deprecated ${var} variable and the new variable SIESTA_${var} are set")
      message(STATUS "  SIESTA_${var}: ${SIESTA_${var}}")
      message(STATUS "  ${var}: ${${var}}")
      message(WARNING " >> The new variable SIESTA_${var} takes precedence")
      set(_dep_var_succeed FALSE)

    endif()
  endif()

endforeach()

list(POP_BACK CMAKE_MESSAGE_INDENT)
if(_dep_var_succeed)
  message(CHECK_PASS "success")
else()
  message(CHECK_FAIL "fail - one or more deprecated variable uses")
endif()

