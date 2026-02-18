# This is the main engine for SIESTA tests.
include(tests/SiestaTests)

list(APPEND CMAKE_MESSAGE_CONTEXT test)
find_package(Python3
  3.5.0 # subprocess limitation
  COMPONENTS Interpreter
  )

add_subdirectory("WaterDimer")
add_subdirectory("LinkAtom-DiBenzene")
add_subdirectory("QMonly-Benzene")
add_subdirectory("MM-Gly")
add_subdirectory("MM-Ala")
add_subdirectory("Restraints")
add_subdirectory("MM-Fixed")
add_subdirectory("Elecfield")

# Pop back the test context
list(POP_BACK CMAKE_MESSAGE_CONTEXT)