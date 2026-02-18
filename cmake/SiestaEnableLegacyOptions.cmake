#------------------------------------------------------
#
#  This file contains interim functionality to support the legacy form of options
#  and other variables (e.g. WITH_MPI instead of the preferred SIESTA_WITH_MPI).
#
# In modern versions of CMake (>= 3.13) the existence of a variable with the same name
# as an option will make the 'option' command a no-op, and the option variable will take the
# pre-existing value, but it will not be put in the cache:
#
#   In CMake 3.13 and above the option() command prefers to do nothing
#   when a normal variable of the given name already exists. It does not
#   create or update a cache entry or remove the normal variable. The
#   new behavior is consistent between the first and later runs in a
#   build tree. This policy provides compatibility with projects that
#   have not been updated to expect the new behavior.
#
# We take advantage of this to implement a simple mechanism to honor old-style
# 'option' names: We copy their values into normal variables with the new name.
#
# Variables set explicitly with the form WITH_XXX will be cache variables,
# and the matching SIESTA_WITH_XXX normal variables.
#
# Variables of the form SIESTA_WITH_YYY set explicitly will be cache variables.
#
# Variables of the form SIESTA_WITH_YYY not set explicitly will also be cache variables, as
# they will be set to default values inside Siesta with the 'option' statement.

# This 'asymmetry' is not bad in principle, and might actually be informative.
# Note also that the cache will in all cases contain the relevant information to reproduce
# all the settings across 'cmake' runs, even if some of the variables are normal, as they
# will be derived in all runs from cache information.
#

message(WARNING
"Allowing non-prefixed variables (e.g. WITH_MPI --> SIESTA_WITH_MPI).\
 This feature might go at any time. Please update your build scripts!")

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
      set(SIESTA_${opt} ${${opt}})
      message(STATUS " >> Using value of deprecated ${opt} option variable for new variable SIESTA_${opt}")
     else ()
      message(STATUS " >> Both the deprecated ${opt} option variable and the new variable SIESTA_${opt} are set")
      message(STATUS "SIESTA_${opt}: ${SIESTA_${opt}}")
      message(WARNING " >> The new variable SIESTA_${opt} takes precedence")
     endif()
   endif()

endforeach()

#       Variables that are not really options.
#
#       SIESTA_WITH_UNIT_CONVENTION is a string with values in a certain enumeration.
#       Note that it is not really checked for allowed values, except if the GUI is used
#       (but it is an "advanced" setting).
#       This variable is set in the top CMakeLists.txt file, and a normal variable will take precedence.
#
#       PROFILE_NVTX_LIBRARY is not really an option, but a path.
#
#       NOTE: WANNIER90_PACKAGE is currently expected to be an environment variable
#
#       If a super-project sets a normal variable of the same name, it will take precedence. It does not
#       need to be made into a cache variable, since the settings will be reproducible across 'cmake' runs.

foreach(var
        WITH_UNIT_CONVENTION
	PROFILE_NVTX_LIBRARY
)

   if( DEFINED ${var} )
     if ( NOT DEFINED SIESTA_${var})
      set(SIESTA_${var} ${${var}})
      message(STATUS " >> Using value of deprecated ${var} variable for new normal variable SIESTA_${var}")
     else ()

      message(STATUS " >> Both the deprecated ${var} variable and the new variable SIESTA_${var} are set")
      message(STATUS "SIESTA_${var}: ${SIESTA_${var}}")
      message(WARNING " >> The new variable SIESTA_${var} takes precedence")
      
     endif()
   endif()

endforeach()

