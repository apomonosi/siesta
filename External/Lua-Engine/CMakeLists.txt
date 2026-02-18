# Define the root folder of flook
set(FLOOK_ROOT "$ENV{FLOOK_ROOT}" CACHE FILEPATH "flook installation path")

# Try first a proper installation with a valid pkg-config file,
# augmenting any specified search paths with FLOOK_ROOT

if(NOT "${FLOOK_ROOT}" STREQUAL "" )
  list(APPEND CMAKE_PREFIX_PATH "${FLOOK_ROOT}")
endif()

find_package(PkgConfig REQUIRED QUIET)

macro(found_return message)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
  set(FLOOK_FOUND TRUE)
  if("${message}" STREQUAL "")
    message(CHECK_PASS "found")
  else()
    message(CHECK_PASS "found ${message}")
  endif()
  return()
endmacro()


message(STATUS "Searching for flook library")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

message(CHECK_START "trying pkg-config")
pkg_check_modules(FLOOK flook>=0.8.4)


macro(flook_add_dependencies)

  # TODO try if dl should be searched for by the libraries
  # this is very hard since we should do all combinations (likely)
  # dl is not a particularly good search name (user-names, etc are likely to
  # hit this key)
  message(CHECK_START "Searching for library [dl]")
  find_library(_flook_lib_dl dl
    DOC "Searching for dl library (dependency on flook)"
    REQUIRED
    )
  message(CHECK_PASS "Found library [dl]: ${_flook_lib_dl}")

  target_link_libraries(flook::flook INTERFACE ${_flook_lib_dl})

  message(CHECK_START "Searching for library [readline]")
  find_library(_flook_lib_readline readline
    DOC "Searching for readline library (dependency on flook)"
    REQUIRED
    )
  target_link_libraries(flook::flook INTERFACE ${_flook_lib_readline})
  message(CHECK_PASS "Found library [readline]: ${_flook_lib_readline}")

  # If the environment variable LUA_DIR is defined, flook will not compile
  # lua, and will not include it in libflookall. Therefore, we need to
  # explicitly include it.
  if(DEFINED ENV{LUA_DIR})
    message(CHECK_START "Searching for library [lua]")
    find_library(_flook_lib_lua lua
      DOC "Searching for lua library (dependency on flook)"
      REQUIRED
      PATHS ENV LUA_DIR
      PATH_SUFFIXES lib lib64
    )
    target_link_libraries(flook::flook INTERFACE ${_flook_lib_lua})
    message(CHECK_PASS "Found library [lua]: ${_flook_lib_lua}")
  endif()
endmacro()


if(FLOOK_FOUND)

  message(VERBOSE "FLOOK_LINK_LIBRARIES: ${FLOOK_LINK_LIBRARIES}")
  message(VERBOSE "FLOOK_INCLUDE_DIRS: ${FLOOK_INCLUDE_DIRS}")

  # Add libraries to namespace
  add_library(flook::flook INTERFACE IMPORTED GLOBAL)

  # add libraries
  target_link_libraries(flook::flook
    INTERFACE "${FLOOK_LINK_LIBRARIES}")
  target_include_directories(flook::flook
    INTERFACE "${FLOOK_INCLUDE_DIRS}")

  found_return("")
else()
  message(CHECK_FAIL "not found")
endif()

# Check for an installation without .pc file in FLOOK_ROOT and set libraries
# manually, assuming the recommended recipe has been followed...
# Search for the library in flook_root
message(CHECK_START "Checking libraries (flookall) in FLOOK_ROOT=${FLOOK_ROOT}")
find_library(FLOOK_LINK_LIBRARIES
  flookall
  HINTS ${FLOOK_ROOT}
  PATH_SUFFIXES lib lib64
  DOC "Libraries for flook"
  )

if("${FLOOK_LINK_LIBRARIES}" STREQUAL "FLOOK_LINK_LIBRARIES-NOTFOUND")
  message(CHECK_FAIL "not found")
  unset(FLOOK_LINK_LIBRARIES CACHE)
else()
  set(FLOOK_INCLUDE_DIRS "${FLOOK_ROOT}/include")
  message(VERBOSE "FLOOK_LINK_LIBRARIES: ${FLOOK_LINK_LIBRARIES}")
  message(VERBOSE "FLOOK_INCLUDE_DIRS: ${FLOOK_INCLUDE_DIRS}")

  add_library(flook::flook INTERFACE IMPORTED GLOBAL)
  # We know we have not added dl library
  target_link_libraries(flook::flook
    INTERFACE "${FLOOK_LINK_LIBRARIES}" "dl")
  target_include_directories(flook::flook
    INTERFACE "${FLOOK_INCLUDE_DIRS}")

  flook_add_dependencies()

  found_return("")
endif()


# Last resort will be to compile flook on the fly

include(ExternalProject)

# Compiling the internal Lua package requires tweaking if you are
# using a different platform than `linux`.
#
#    PLATFORM = aix | bsd | c89 | freebsd | generic
#               | linux | macosx | mingw | posix | solaris
#
# where the makefile system will try and guess the correct
# platform. However, if the build fails due to the Lua library, then
# please supply `PLATFORM` with the correct platform in your
# `setup.make` file.
#
# Please add more cases here if needed.

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
  set(FLOOK_PLATFORM "macosx")
else()
  set(FLOOK_PLATFORM "linux")
endif()

configure_file(
  setup.make.in
  "${CMAKE_CURRENT_BINARY_DIR}/setup.make"
  @ONLY
  )

# Find the make-command
find_program(MAKE_EXECUTABLE NAMES gmake make mingw32-make REQUIRED)


# Possible sources
set(url_default "https://github.com/ElectronicStructureLibrary/flook/archive/refs/tags/v0.8.4.tar.gz")
#
if(NOT "$ENV{FLOOK_PACKAGE}" STREQUAL "")

  set(url "$ENV{FLOOK_PACKAGE}")
  message(STATUS "flook: will compile on-the-fly from sources: $ENV{FLOOK_PACKAGE}")

elseif(EXISTS "${PROJECT_SOURCE_DIR}/External/Lua-Engine/flook/LICENSE")
  
  # The sub-module exists, or at least the LICENSE file does!

  set(url "${PROJECT_SOURCE_DIR}/External/Lua-Engine/flook")
  set(FLOOK_SOURCE_DIR "${url}" PARENT_SCOPE)
  message(STATUS "flook: will compile on-the-fly from submodule in External/Lua-Engine/flook")

else()

    set(url "${url_default}")

    message(WARNING
      "flook: Trying to download the sources directly from:"
      "  ${url_default} "
      "If this is not possible (no access to internet) then please define the "
      "environment variable FLOOK_PACKAGE pointing to a flook tarball with the sources. "
      "The easiest solution would be to check-out the submodule with the following commands:\n"
      "  cd <siesta-top-directory>\n"
      "  git submodule update --init External/Lua-Engine/flook"
      )
endif()

set(_flook_build_dir
  ${CMAKE_CURRENT_BINARY_DIR}/flook-prefix/src/flook-build)

set(FLOOK_INCLUDE_DIRS
  ${_flook_build_dir})
set(FLOOK_LINK_LIBRARIES
  ${_flook_build_dir}/libflookall.a)

ExternalProject_Add(flook
 URL "${url}"
 UPDATE_DISCONNECTED true
 ##CONFIGURE_HANDLED_BY_BUILD true
 DOWNLOAD_EXTRACT_TIMESTAMP TRUE
 CONFIGURE_COMMAND ""
 BUILD_COMMAND ${MAKE_EXECUTABLE} -j liball
 ## We avoid installing until we know how to specify the prefix...
 ##  ... the following should work but does not.
 ## INSTALL_COMMAND ${MAKE_EXECUTABLE} DESTDIR=<INSTALL_DIR> install
 INSTALL_COMMAND ""
 TEST_COMMAND ""
 #
 # See https://stackoverflow.com/questions/54866067/cmake-and-ninja-missing-and-no-known-rule-to-make-it
 #
 BUILD_BYPRODUCTS ${FLOOK_LINK_LIBRARIES}
)

#
# Copy configured setup.make file
# Copy also everything from SOURCE_DIR to BINARY_DIR (this is
# typically done by 'configure' for autotools, but it is not
# done for 'make' projects -- there must be a better way).
# Note special idiom for the 'cp' command. **Check** that it
# works with all shells?
#
ExternalProject_Add_Step(flook
  copy_setup_make
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/setup.make . && cp -a <SOURCE_DIR>/. .
  DEPENDEES download
  DEPENDERS build
  WORKING_DIRECTORY <BINARY_DIR>
)

add_library(flook::flook INTERFACE IMPORTED GLOBAL)

add_dependencies(flook::flook flook)

target_link_libraries(flook::flook
  INTERFACE
  ${FLOOK_LINK_LIBRARIES}
  )

target_include_directories(flook::flook
  INTERFACE
  ${FLOOK_INCLUDE_DIRS}
)

flook_add_dependencies()
