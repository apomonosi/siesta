# find CheSS via brute force for now
# (wip: Not tested: some tweaks might be still needed)
#
# The environment variable CHESS_ROOT must point to the
# root of the CheSS installation

message(CHECK_START "Searching for CheSS library")

include(FindPackageHandleStandardArgs)

# The (static) libraries are all in ${CHESS_ROOT}/lib
# ... but we provide extra hooks for other places
# for the subordinate libraries
#
# WARNING: Be sure to check that some things
# are not picked up from other places, such as a PSOLVER installation...
# The following deactivates searching in the default places

set(avoided_heuristics
    NO_CMAKE_PATH
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    )

find_library(CHESS_LIBRARIES
  ${avoided_heuristics}
  NAMES CheSS-1
  PATH_SUFFIXES lib lib64
  HINTS
  ENV CHESS_ROOT
  DOC "CheSS libraries list")

find_library(FUTILE_CHESS_LIBRARIES
  ${avoided_heuristics}
  NAMES futile-1
  PATH_SUFFIXES lib lib64
  HINTS
  ENV CHESS_ROOT
  ENV FUTILE_CHESS_ROOT
  DOC "Futile library for use with CheSS")

find_library(YAML_CHESS_LIBRARIES
  ${avoided_heuristics}
  NAMES yaml
  PATH_SUFFIXES lib lib64
  HINTS
  ENV CHESS_ROOT
  ENV YAML_CHESS_ROOT
  DOC "yaml library for use with CheSS")

find_path(CHESS_INCLUDE_DIR
  ${avoided_heuristics}
  NAMES foe_base.mod   # ** CHECK
  PATH_SUFFIXES include
  HINTS
  ENV CHESS_ROOT)

find_package_handle_standard_args(CustomCHESS "DEFAULT_MSG"
                                  CHESS_LIBRARIES FUTILE_CHESS_LIBRARIES 
				  YAML_CHESS_LIBRARIES CHESS_INCLUDE_DIR)


if(NOT CustomCHESS_FOUND)
  message(CHECK_FAIL "not found")
  return()
endif()

message("CHESS_INCLUDE_DIR: ${CHESS_INCLUDE_DIR}")
message("CHESS_LIBRARIES: ${CHESS_LIBRARIES}")
message("FUTILE_CHESS_LIBRARIES: ${FUTILE_CHESS_LIBRARIES}")
message("YAML_CHESS_LIBRARIES: ${YAML_CHESS_DIST_LIBRARIES}")


if(CustomCHESS_FOUND AND NOT TARGET CheSS::CheSS)


# Instead of 'STATIC' 
# we might want to use 'SHARED' if appropriate (to be implemented)

  add_library(chess STATIC IMPORTED)
  set_target_properties(chess PROPERTIES
  			      IMPORTED_LOCATION "${CHESS_LIBRARIES}")
			      
  add_library(futile_chess STATIC IMPORTED)
  set_target_properties(futile_chess PROPERTIES
  			      IMPORTED_LOCATION "${FUTILE_CHESS_LIBRARIES}")

  add_library(yaml_chess STATIC IMPORTED)
  set_target_properties(yaml_chess PROPERTIES
  			     IMPORTED_LOCATION "${YAML_CHESS_DIST_LIBRARIES}")

  add_library(CheSS::CheSS INTERFACE IMPORTED)

  target_link_libraries(CheSS::CheSS
                          INTERFACE
                          chess
			  futile_chess
			  yaml_chess)
			  
  set_target_properties(CheSS::CheSS
                        PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${CHESS_INCLUDE_DIR}")

endif()

message(CHECK_PASS "found")

