include(SiestaFindPackage)

Siesta_find_package(mstore
  REQUIRED
  GIT_REPOSITORY "https://github.com/grimme-lab/mstore.git"
  GIT_TAG "HEAD"
  SOURCE_DIR "${DFTD3_SOURCE_ROOT}/mstore"
  )
