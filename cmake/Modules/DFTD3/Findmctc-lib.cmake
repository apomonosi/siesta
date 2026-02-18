include(SiestaFindPackage)

Siesta_find_package(mctc-lib
  REQUIRED
  GIT_REPOSITORY "https://github.com/grimme-lab/mctc-lib.git"
  GIT_TAG "v0.3.1"
  SOURCE_DIR "${DFTD3_SOURCE_ROOT}/mctc-lib"
  )
