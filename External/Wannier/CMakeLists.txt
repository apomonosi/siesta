if (NOT TARGET wrapper-wannier90::wrapper-wannier90)

  set(CMAKE_PROJECT_wrapper-wannier90_INCLUDE_BEFORE
         "${PROJECT_SOURCE_DIR}/External/Options/options_wrapper-wannier90.cmake")
	 
  # Location of wannier90-3.1.0.tar.gz package
  # It can be a local file or a remote URL
  
  if(NOT "$ENV{WANNIER90_PACKAGE}" STREQUAL "")
       set(url "$ENV{WANNIER90_PACKAGE}")
  elseif(EXISTS
     "${PROJECT_SOURCE_DIR}/External/Wannier/wrapper-wannier90/LICENSE")
     #
     # We assume that we have somehow put the pristine source in this directory
     # (perhaps using a submodule -- but do not 'init' it by default,
     # as it is huge -- cloned directory is over 170Mb)
     # Better ship tar.gz file around
     #
     set(url "${PROJECT_SOURCE_DIR}/External/Wannier/wrapper-wannier90")

  elseif(EXISTS
     "${PROJECT_SOURCE_DIR}/External/Wannier/patched_wrapper-wannier90/CMakeLists.txt")

    # Intended for developers
    message(STATUS "Adding External/Wannier/patched_wrapper-wannier90 subdirectory...")
    message(STATUS "... containing already-patched sources for w90 wrapper")

    add_subdirectory("${PROJECT_SOURCE_DIR}/External/Wannier/patched_wrapper-wannier90")
    add_library("wrapper-wannier90::wrapper-wannier90" ALIAS "wrapper-wannier90-lib")
    RETURN()
    
  else()

    set(url "https://github.com/wannier-developers/wannier90/archive/v3.1.0.tar.gz")
    
    message(STATUS "About to download the Wannier90 source package")
    message(STATUS "If download is not possible, or too heavy:")
    message(STATUS " Please set environment variable WANNIER90_PACKAGE")
    message(STATUS " to point to a pristine (possibly remote) wannier90-3.1.0.tar.gz package,")
    message(STATUS " or populate ${PROJECT_SOURCE_DIR}/External/Wannier/wrapper-wannier90/")
    message(STATUS " with a pristine wannier90 v3.1.0 src/ directory plus the LICENSE file,")
    message(STATUS " or populate ${PROJECT_SOURCE_DIR}/External/Wannier/patched_wrapper-wannier90/")
    message(STATUS " with the already patched w90 wrapper sources")

  endif()
  
  set(package "wrapper-wannier90")
  message(STATUS "Retrieving pristine sources for wannier90-3.1.0 from ${url}")
  include(FetchContent)

  # Common declaration for both approaches
  FetchContent_Declare(wrapper-wannier90
    DOWNLOAD_EXTRACT_TIMESTAMP true
    PATCH_COMMAND patch -p1 -i ${PROJECT_SOURCE_DIR}/External/Wannier/Patches/3.1.0.patch
    UPDATE_DISCONNECTED true
    UPDATE_COMMAND ""
    URL "${url}"
    URL_HASH SHA256=40651a9832eb93dec20a8360dd535262c261c34e13c41b6755fa6915c936b254
  )

  # Use modern approach for CMake 3.24+ (when FetchContent_MakeAvailable got more robust)
  # Use traditional approach for older versions to ensure patch safety
  if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
    message(STATUS "Using modern FetchContent_MakeAvailable approach (CMake ${CMAKE_VERSION})")
    FetchContent_MakeAvailable(wrapper-wannier90)
  else()
    # Slightly indirect way to avoid re-patching if CMake is invoked again.
    # The call to _Populate, as opposed to _MakeAvailable, seems not to
    # apply the patching step a second time...
    # (set -DFETCHCONTENT_QUIET=OFF to see the internal messages)
    message(STATUS "Using traditional FetchContent approach for compatibility (CMake ${CMAKE_VERSION})")
    
    FetchContent_GetProperties(wrapper-wannier90)
    if(NOT wrapper-wannier90_POPULATED)
      # Fetch the content using previously declared details
      message(STATUS "Making-available ${package}...")
      FetchContent_Populate(wrapper-wannier90)

      # Set custom variables, policies, etc.
      # ...
    endif()

    # Bring the populated content into the build
    message(STATUS "Adding subdirectory ${wrapper-wannier90_SOURCE_DIR}")
    add_subdirectory("${wrapper-wannier90_SOURCE_DIR}" "${wrapper-wannier90_BINARY_DIR}")
  endif()

  # The wrapper-wannier90::wrapper-wannier90 target is not defined in the
  # subproject until installation time. Hence the need for the
  # following lines.
  # Note the GLOBAL qualifier, since this is in a subordinate scope
  # I am not sure about the interpretation of 'wrapper-wannier90'
  # in the second line.
  
  add_library(wrapper-wannier90::wrapper-wannier90 INTERFACE IMPORTED GLOBAL)
  target_link_libraries(wrapper-wannier90::wrapper-wannier90 INTERFACE wrapper-wannier90)

  # We need the module directory in the subproject before we finish the configure stage
  FetchContent_GetProperties(${package} BINARY_DIR WRAPPER_WANNIER90_BINARY_DIR)
  if(NOT EXISTS "${WRAPPER_WANNIER90_BINARY_DIR}/include")
     make_directory("${WRAPPER_WANNIER90_BINARY_DIR}/include")
  endif()

endif()
