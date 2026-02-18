include(tests/SiestaTests)

set(LINDHARD_INP_DIR "${SIESTA_UTILS_SOURCE_DIR}/Lindhard/Test")
set(LINDHARD_REF_DIR "${SIESTA_UTILS_SOURCE_DIR}/Lindhard/Test/Reference")
set(LINDHARD_BIN_DIR "${SIESTA_UTILS_BINARY_DIR}/Lindhard/Src")

#Utilities
siesta_util_test(
  NAME       "lindhard"
  BASH_SCR   "${LINDHARD_INP_DIR}/execute_test.sh"
  BASH_VERIF "${LINDHARD_INP_DIR}/verify_test.sh"
  EXECUTABLE "lindhard"
  INPUT_FILE "ldau.fdf"
  EXEC_DIR   "${LINDHARD_BIN_DIR}"
  SOURCE_DIR "${LINDHARD_INP_DIR}"
  OUTFILE    "lindhard.out"
  REFDIR     "${LINDHARD_REF_DIR}"
)