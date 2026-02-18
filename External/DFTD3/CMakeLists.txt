#
# s-dftd3 depends on a number of libraries
# They can be included as submodules in this directory,
# and this file takes care of the details.
#
# Originally we processed the dependencies in order, but
# now we just look for s-dftd3 and let its own build system
# (or cmake config file) start a search for the dependencies.
#
# It would be good to have an option to build only
# the library itself.

# Configuration of sub-project variables ----------------

# (Will be set to "" if the SIESTA option is not defined, thus
#  turning the 'option' setting inside DFTD3 into a no-op,
#  and effectively disabling the tests.)
# (If quotes are not used, the behavior might change.)
# set(WITH_TESTS "${SIESTA_WITH_DFTD3_TESTS}")
#
# It is clearer to put this here to set a default
#
option(SIESTA_WITH_DFTD3_TESTS "DFT-D3 tests" OFF)
set(WITH_TESTS ${SIESTA_WITH_DFTD3_TESTS})

#
# (How do we turn off the auxiliary tests for toml-f et al?)
#

# We need to forcefully disable OpenMP, s-dftd3 is defaulting to true for this.
set(WITH_OpenMP "${SIESTA_WITH_OPENMP}")

# -------------------------------------------------------

### Re-use information already found by Siesta itself
### and skip BLAS search in s-dftd3

if (TARGET BLAS::BLAS)
  # Re-use target
  set(BLAS_FOUND "TRUE")
elseif(TARGET LAPACK::LAPACK)
  # Create a BLAS wrapper target
  add_library("BLAS::BLAS" INTERFACE IMPORTED)
  target_link_libraries("BLAS::BLAS" INTERFACE LAPACK::LAPACK)
  set(BLAS_FOUND "TRUE")
endif()


# Just search for s-dftd3 ----------------------------------
#
# Use our own finders for possible s-dftd3 dependencies, as we can have source directories
# in places different from the expectation of s-dftd3's finders.
#
list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/DFTD3")
message(DEBUG "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

# Define the source directory root
set(DFTD3_SOURCE_ROOT "${PROJECT_SOURCE_DIR}/External/DFTD3")

if (NOT TARGET s-dftd3::s-dftd3)

  # The REQUIRED flag might not be passed to the module... check
  find_package(Custom-s-dftd3 REQUIRED MODULE)

  if (TARGET s-dftd3::s-dftd3)
    # This is necessary in case a cmake package is found
    # for s-dftd3, as imported targets are not global by
    # default, and we are running this in a "subdirectory"
    # scope. If s-dftd3 is compiled ("source" method), this would not be
    # necessary, since we manually add the library as IMPORTED GLOBAL
    # in SiestaFindPackage.
    #
    set_target_properties(s-dftd3::s-dftd3
                          PROPERTIES
                          IMPORTED_GLOBAL TRUE)
  endif()
endif()
list(POP_FRONT CMAKE_MODULE_PATH)
# ---------------------------------------

