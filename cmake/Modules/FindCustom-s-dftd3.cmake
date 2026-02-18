include(SiestaFindPackage)

Siesta_find_package(s-dftd3
  REQUIRED
  GIT_REPOSITORY "https://github.com/dftd3/simple-dftd3.git"
  GIT_TAG "v1.1.0"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/External/DFTD3/s-dftd3"
  )
