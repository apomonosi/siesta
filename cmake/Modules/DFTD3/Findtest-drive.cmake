include(SiestaFindPackage)

Siesta_find_package(test-drive
  REQUIRED
  GIT_REPOSITORY "https://github.com/fortran-lang/test-drive.git"
  GIT_TAG "v0.4.0"
  SOURCE_DIR "${DFTD3_SOURCE_ROOT}/test-drive"
  )
