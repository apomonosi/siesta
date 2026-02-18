#APostnikov Utilities to test
include(tests/SiestaTests)

set(APOSTNIKOV_INP_DIR "${SIESTA_UTILS_SOURCE_DIR}/Contrib/APostnikov/Tests/Inputs")
set(APOSTNIKOV_REF_DIR "${SIESTA_UTILS_SOURCE_DIR}/Contrib/APostnikov/Tests/Reference")
set(APOSTNIKOV_BIN_DIR "${SIESTA_UTILS_BINARY_DIR}/Contrib/APostnikov")

#Utilities
siesta_util_test(
  NAME       "eig2bxsf"
  BASH_SCR   "${APOSTNIKOV_INP_DIR}/execute_test.sh"
  BASH_VERIF "${APOSTNIKOV_INP_DIR}/verify_test.sh"
  EXECUTABLE "eig2bxsf"
  INPUT_FILE "${APOSTNIKOV_INP_DIR}/eig2bxsf_input.txt"
  EXEC_DIR   "${APOSTNIKOV_BIN_DIR}"
  SOURCE_DIR "${APOSTNIKOV_INP_DIR}"
  OUTFILE    "eig2bxsf.LOG"
  REFDIR     "${APOSTNIKOV_REF_DIR}"
)

siesta_util_test(
  NAME       "fmpdos"
  BASH_SCR   "${APOSTNIKOV_INP_DIR}/execute_test.sh"
  BASH_VERIF "${APOSTNIKOV_INP_DIR}/verify_test.sh"
  EXECUTABLE "fmpdos"
  INPUT_FILE "${APOSTNIKOV_INP_DIR}/fmpdos_input.txt"
  EXEC_DIR   "${APOSTNIKOV_BIN_DIR}"
  SOURCE_DIR "${APOSTNIKOV_INP_DIR}"
  OUTFILE    "diamond.PDOS.fmpdos"
  REFDIR     "${APOSTNIKOV_REF_DIR}"
)

siesta_util_test(
  NAME       "md2axsf"
  BASH_SCR   "${APOSTNIKOV_INP_DIR}/execute_test.sh"
  BASH_VERIF "${APOSTNIKOV_INP_DIR}/verify_test.sh"
  EXECUTABLE "md2axsf"
  INPUT_FILE "${APOSTNIKOV_INP_DIR}/md2axsf_input.txt"
  EXEC_DIR   "${APOSTNIKOV_BIN_DIR}"
  SOURCE_DIR "${APOSTNIKOV_INP_DIR}"
  OUTFILE    "diamond.AXSF"
  REFDIR     "${APOSTNIKOV_REF_DIR}"
)

siesta_util_test(
  NAME       "xv2vesta"
  BASH_SCR   "${APOSTNIKOV_INP_DIR}/execute_test.sh"
  BASH_VERIF "${APOSTNIKOV_INP_DIR}/verify_test.sh"
  EXECUTABLE "xv2vesta"
  INPUT_FILE "${APOSTNIKOV_INP_DIR}/xv2vesta_input.txt"
  EXEC_DIR   "${APOSTNIKOV_BIN_DIR}"
  SOURCE_DIR "${APOSTNIKOV_INP_DIR}"
  OUTFILE    "diamond.vesta"
  REFDIR     "${APOSTNIKOV_REF_DIR}"
)

# the rho2xsf util has to be before the xv2xsf because both utils outputs the
# same file.
siesta_util_test(
  NAME       "rho2xsf"
  BASH_SCR   "${APOSTNIKOV_INP_DIR}/execute_test.sh"
  BASH_VERIF "${APOSTNIKOV_INP_DIR}/verify_test.sh"
  EXECUTABLE "rho2xsf"
  INPUT_FILE "${APOSTNIKOV_INP_DIR}/rho2xsf_input.txt"
  EXEC_DIR   "${APOSTNIKOV_BIN_DIR}"
  SOURCE_DIR "${APOSTNIKOV_INP_DIR}"
  OUTFILE    "diamond.XSF"
  REFDIR     "${APOSTNIKOV_REF_DIR}"
)

# This test has to be after the rho2xsf because xv2xsf can rewrite the output
# without problems. Both Utils output the same file name.
siesta_util_test(
  NAME       "xv2xsf"
  BASH_SCR   "${APOSTNIKOV_INP_DIR}/execute_test.sh"
  BASH_VERIF "${APOSTNIKOV_INP_DIR}/verify_test.sh"
  EXECUTABLE "xv2xsf"
  INPUT_FILE "${APOSTNIKOV_INP_DIR}/xv2xsf_input.txt"
  EXEC_DIR   "${APOSTNIKOV_BIN_DIR}"
  SOURCE_DIR "${APOSTNIKOV_INP_DIR}"
  OUTFILE    "diamond.XSF"
  REFDIR     "${APOSTNIKOV_REF_DIR}"
)