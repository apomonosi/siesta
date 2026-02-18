# utility package to detect conda builds

set(SIESTA_conda_build "$ENV{CONDA_BLD_PATH}")


if(NOT SIESTA_conda_build)
  message(STATUS "Normal build, no guard against relocation of prefix")

  macro(siesta_conda_fix var)
  endmacro()

  return()
endif()

message(STATUS "Building in conda CI, replacing $PREFIX,$BUILD_PREFIX,$SRC_DIR")

macro(siesta_conda_fix var)
  string(REPLACE "$PREFIX"
    "<prefix>"
    ${var} "${${var}}"
    )
  string(REPLACE "$ENV{PREFIX}"
    "<prefix>"
    ${var} "${${var}}"
    )
  string(REPLACE "$BUILD_PREFIX"
    "<prefix>/_build_env"
    ${var} "${${var}}"
    )
  string(REPLACE "$ENV{BUILD_PREFIX}"
    "<prefix>/_build_env"
    ${var} "${${var}}"
    )
  string(REPLACE "$SRC_DIR"
    "<prefix>/work"
    ${var} "${${var}}"
    )
  string(REPLACE "$ENV{SRC_DIR}"
    "<prefix>/work"
    ${var} "${${var}}"
    )
endmacro()
