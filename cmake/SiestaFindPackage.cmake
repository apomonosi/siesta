# This file is a modification of the one with the
# name "sdftd3_find_package" in s-dft3. The original terms follow:
# -----------------------------------------
# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

# This changed version streamlines the effects and uses a function rather
# than a macro.
# It enables a bit more things and ensures expressivity in the options.
# Upon return it returns in ${NAME}_FOUND_METHOD the method by which it
# was found (only if ${NAME}_FOUND is true).

# Handling of subproject dependencies
# This function should be used in FindXXX.cmake files since
# it directly utilizes the _FIND_REQUIRED|QUIETLY variables
# This is important to not duplicate efforts.
function(Siesta_find_package)
  # Define options
  set(options
    REQUIRED)
  set(oneValueArgs
    # the package name:
    # Used for variable lookups {NAME}_FIND_METHODS etc.
    NAME
    MIN_VERSION
    GIT_REPOSITORY GIT_TAG # for git repo's
    TARGET # name of the created target (defaluts to NAME::NAME)
    SOURCE_DIR)
  set(multiValueArgs
    CMAKE # name(s) of cmake package (default to NAME)
    PKG_CONFIG # name(s) of pkg-config package (default to NAME)
    )
  cmake_parse_arguments(_f "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # check arguments
  if(NOT DEFINED _f_NAME)
    # Check that the package as an unnamed argument
    list(LENGTH _f_UNPARSED_ARGUMENTS nargs)
    if(nargs GREATER 0)
      list(POP_FRONT _f_UNPARSED_ARGUMENTS _f_NAME)
    else()
      message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION} missing argument: NAME")
    endif()
  endif()

  # Define the CMAKE and PKG_CONFIG variables
  if(NOT DEFINED _f_CMAKE)
    set(_f_CMAKE "${_f_NAME}")
  endif()
  if(NOT DEFINED _f_PKG_CONFIG)
    set(_f_PKG_CONFIG "${_f_NAME}")
  endif()
  if(NOT DEFINED _f_TARGET)
    set(_f_TARGET "${_f_NAME}::${_f_NAME}")
  endif()


  set(pkg "${_f_NAME}")
  # convert package name to lower and upper case
  string(TOLOWER "${pkg}" pkg_lc)
  string(TOUPPER "${pkg}" pkg_uc)

  # Default to not have been found.
  # These are looked up in the loop constructs, so we
  # will change them while running (and in the end)
  set("${pkg}_FOUND" FALSE PARENT_SCOPE)
  set("${pkg_uc}_FOUND" FALSE PARENT_SCOPE)

  macro(mymsg check msg)
    if(NOT ${_f_NAME}_FIND_QUIETLY)
      message(${check} "${msg}")
    endif()
  endmacro()

  mymsg(CHECK_START "Searching for ${pkg}")
  # always add space (for sub-stuff)
  list(APPEND CMAKE_MESSAGE_INDENT "| ")

  # Allow users to pass options via flags
  #  *_GIT_REPOSITORY
  #  *_GIT_TAG
  #  *_SOURCE_DIR
  set(_f_REPO "${_f_GIT_REPOSITORY}")
  set(_f_TAG "${_f_GIT_TAG}")
  foreach(n IN ITEMS ${pkg_lc} ${pkg_uc} ${pkg})

    # Also allow externally setting the GIT_REPOSITORY and GIT_TAG
    if(DEFINED "${n}_GIT_REPOSITORY")
      set(_f_REPO "${${n}_GIT_REPOSITORY}")
    endif()

    if(DEFINED "${n}_GIT_TAG")
      set(_f_TAG "${${n}_GIT_TAG}")
    endif()

    # The submodule directory is a SOURCE_DIR from the outside
    # world. Hence it should look up the SOURCE_DIR
    if(DEFINED "${n}_SOURCE_DIR")
      set(_f_SOURCE_DIR "${${n}_SOURCE_DIR}")
    endif()

    mark_as_advanced(${n}_GIT_REPOSITORY ${n}_GIT_TAG ${n}_SOURCE_DIR)
  endforeach()


  # Determine the allowed find-methods
  set(all_f_methods "cmake" "pkgconf" "source" "fetch")
  set(allowed_f_methods "cmake" "pkgconf")
  if(DEFINED _f_SOURCE_DIR)
    list(APPEND allowed_f_methods "source")
  endif()
  if(DEFINED "_f_REPO" AND DEFINED "_f_TAG")
    list(APPEND allowed_f_methods "fetch")
  elseif(DEFINED "_f_REPO")
    message(WARNING "${CMAKE_CURRENT_FUNCTION}: ${pkg} will not allow fetch as the GIT_REPOSITORY is supplied, but no GIT_TAG, please use -D${pkg}_GIT_TAG=")
  elseif(DEFINED "_f_TAG")
    message(WARNING "${CMAKE_CURRENT_FUNCTION}: ${pkg} will not allow fetch as the GIT_TAG is supplied, but no GIT_REPOSITORY, please use -D${pkg}_GIT_REPOSITORY=")
  endif()

  # default the searched methods
  set(f_methods "${allowed_f_methods}")
  if(DEFINED "${PROJECT_NAME}_FIND_METHOD")
    set(f_methods "${${PROJECT_NAME}_FIND_METHOD}")
  endif()

  foreach(n IN ITEMS ${pkg_lc} ${pkg_uc} ${pkg})
    # Parse the option of the *find-methods*
    if(DEFINED ${n}_FIND_METHOD)
      set(f_methods "${${n}_FIND_METHOD}")
    endif()
    mark_as_advanced(${n}_FIND_METHOD)

    # Also check if the FOUND method is defined from the cache
    # This should overwrite everything
    if(DEFINED ${n}_FOUND_METHOD)
      set(f_methods "${${n}_FOUND_METHOD}")
      break()
    endif()
  endforeach()

  # Strip away elements that are not in the allowed-f-methods
  # If there is no submodule, it shouldn't be allowed.
  list(REMOVE_DUPLICATES f_methods)
  # lower
  string(TOLOWER "${f_methods}" f_methods)
  foreach(method IN LISTS f_methods)
    if(NOT method IN_LIST allowed_f_methods)
      list(REMOVE_ITEM f_methods ${method})
    endif()
  endforeach()

  if(TARGET "${_f_TARGET}")
    # skip all methods
    message(DEBUG "${CMAKE_CURRENT_FUNCTION}: target already found, no searching done.")
    set(f_methods "")
  endif()


  # A placeholder for what we actually used to find the package
  set(f_method)

  message(STATUS "${CMAKE_CURRENT_FUNCTION}[${pkg}] METHODS | ALLOWED = ${f_methods} | ${allowed_f_methods}")

  foreach(method IN ITEMS ${f_methods})
    # assert that method is in all_f_methods
    if(NOT method IN_LIST all_f_methods)
      message(FATAL_ERROR "${CMAKE_CURRENT_FUNCTION}: find-method ${method} is not in allowed methods ${all_f_methods}")
    endif()

    if("${method}" STREQUAL "cmake")

      foreach(_cmake IN LISTS _f_CMAKE)

        mymsg(CHECK_START "CMake package lookup [${_cmake}]")

        if (DEFINED _f_MIN_VERSION)
          set(_version_string ${_f_MIN_VERSION})
        endif()

        find_package("${_cmake}" ${_version_string} CONFIG QUIET)

        if( ${_cmake}_FOUND OR TARGET "${_f_TARGET}")
          # signal to outer loop
          set(${pkg}_FOUND TRUE)

          # try and figure out targets
          if(NOT TARGET "${_f_TARGET}")

            # we have to point to the other things
            if(TARGET "${_cmake}::${_cmake}")
              add_library("${_f_TARGET}" INTERFACE IMPORTED GLOBAL)
              target_link_libraries("${_f_TARGET}" INTERFACE "${_cmake}::${_cmake}")

            elseif(TARGET "${_cmake}")
              add_library("${_f_TARGET}" INTERFACE IMPORTED GLOBAL)
              target_link_libraries("${_f_TARGET}" INTERFACE "${_cmake}")

            else()
              message(STATUS "Located [${pkg}] through CMake package=${_cmake} but could not figure out target")
              set(${pkg}_FOUND FALSE)
            endif()

          endif()

          if(${pkg}_FOUND)
            mymsg(CHECK_PASS "found")

            # Propagate the _DIR for the package found, so that it can be recorded
            # in the cmake config file

            set(${pkg}_DIR "${${pkg}_DIR}" PARENT_SCOPE)
            message(STATUS "Storing ${pkg}_DIR: ${${pkg}_DIR}")

            break()
          endif()

        endif()

        mymsg(CHECK_FAIL "not found/compatible")
      endforeach()

      if(${pkg}_FOUND)
        set(f_method "cmake")
        break()
      endif()

    endif()

    if("${method}" STREQUAL "pkgconf")
      find_package(PkgConfig QUIET REQUIRED)


      foreach(pkg_config IN LISTS _f_PKG_CONFIG)

        # extract the pkg-config details
        # format:
        #  NAME1 pack1 pack2;NAME2 pack3 pack4 ...
        separate_arguments(pkg_config)
        list(POP_FRONT pkg_config pkg_name)
        list(LENGTH pkg_config pkg_config_len)
        if(${pkg_config_len} EQUAL 0)
          # same package name as the package itself
          set(pkg_config ${pkg_name})
        endif()

        mymsg(CHECK_START "pkg-config package lookup[${pkg_name}]")

        if (DEFINED _f_MIN_VERSION)
	   set(_version_string ">=${_f_MIN_VERSION}")
        endif()

        pkg_check_modules(${pkg_name}
          IMPORTED_TARGET GLOBAL
          ${pkg_config}${_version_string}
          )

        if( ${pkg_name}_FOUND OR TARGET PkgConfig::${pkg_name})
          mymsg(CHECK_PASS "found")


          #add_library("${_f_TARGET}" ALIAS PkgConfig::${pkg_name})
          add_library("${_f_TARGET}" INTERFACE IMPORTED GLOBAL)
          target_link_libraries("${_f_TARGET}" INTERFACE "PkgConfig::${pkg_name}")

          set(${pkg}_FOUND TRUE)

          # Propagate the PREFIX for the package found, so that it can be recorded
          # in the cmake config file

          get_filename_component(${pkg_name}_PREFIX "${${pkg_name}_LIBRARY_DIRS}" DIRECTORY)
          message(STATUS "${pkg_name}_PREFIX: ${${pkg_name}_PREFIX}")
          set(${pkg_name}_PREFIX "${${pkg_name}_PREFIX}" PARENT_SCOPE)

          break()
        else()
          mymsg(CHECK_FAIL "not found/compatible")
        endif()

      endforeach()

      if(${pkg}_FOUND)
        set(f_method "pkgconf")
        break()
      endif()

    endif()

    if("${method}" STREQUAL "source")
      mymsg(CHECK_START "source in folder: ${_f_SOURCE_DIR}")

      set("${pkg_uc}_SOURCE_DIR" "${_f_SOURCE_DIR}")
      if(EXISTS "${${pkg_uc}_SOURCE_DIR}/CMakeLists.txt")
        mymsg(CHECK_PASS "found")

        add_subdirectory(
          "${${pkg_uc}_SOURCE_DIR}"
          "${PROJECT_BINARY_DIR}/${pkg}"
        )

        # Needed in some cases (see below)
	if(NOT EXISTS "${PROJECT_BINARY_DIR}/${pkg}/include")
          make_directory("${PROJECT_BINARY_DIR}/${pkg}/include")
        endif()

        add_library("${_f_TARGET}" INTERFACE IMPORTED GLOBAL)
        target_link_libraries("${_f_TARGET}" INTERFACE ${pkg})

        set(f_method "source")
        break()
      else()
        mymsg(CHECK_FAIL "not found/compatible")
      endif()
    endif()

    if("${method}" STREQUAL "fetch")
      mymsg(CHECK_START "fetching from ${_f_REPO}")

      include(FetchContent)
      FetchContent_Declare("${pkg}"
        GIT_REPOSITORY "${_f_REPO}"
        GIT_TAG "${_f_TAG}"
      )
      FetchContent_MakeAvailable("${pkg}")

      add_library("${_f_TARGET}" INTERFACE IMPORTED GLOBAL)
      target_link_libraries("${_f_TARGET}" INTERFACE "${pkg}")

      # We need the module directory in the subproject before we finish the configure stage
      # Note that this is relevant if the "Fortran_MODULE_DIRECTORY" property is set to
      # the 'include' subdirectory, as it is done for some libraries.
      #
      # (This is found experimentally; it does not seem to apply always, but it does not
      #  hurt to do it.)

      FetchContent_GetProperties("${pkg}" BINARY_DIR "${pkg_uc}_BINARY_DIR")
      message(STATUS "BINARY_DIR for fetched ${pkg}: ${${pkg_uc}_BINARY_DIR}")

      if(NOT EXISTS "${${pkg_uc}_BINARY_DIR}/include")
        make_directory("${${pkg_uc}_BINARY_DIR}/include")
      endif()

      mymsg(CHECK_PASS "fetched")
      set(f_method "fetch")

      break()
    endif()

  endforeach()

  list(POP_BACK CMAKE_MESSAGE_INDENT)
  if(TARGET "${_f_TARGET}")
    mymsg(CHECK_PASS "found")

    foreach(n IN ITEMS ${pkg} ${pkg_lc} ${pkg_uc})
      # notify about the found-variables
      set("${n}_FOUND" TRUE PARENT_SCOPE)
      set("${n}_FOUND_METHOD" ${f_method} CACHE STRING "How the package ${pkg} was found")
    endforeach()

  else()
    mymsg(CHECK_FAIL "not found")

    foreach(n IN ITEMS ${pkg} ${pkg_lc} ${pkg_uc})
      # notify about the found-variables
      set("${n}_FOUND" FALSE PARENT_SCOPE)
    endforeach()

    if( ${_f_NAME}_FIND_REQUIRED OR _f_REQUIRED )
      # package hasn't been found but REQUIRED. Emit an error.
      message(FATAL_ERROR "Required package ${_f_NAME} cannot be found")
    endif()

  endif()

endfunction()

