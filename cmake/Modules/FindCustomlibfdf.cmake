include(SiestaFindPackage)

Siesta_find_package(libfdf
  REQUIRED
  MIN_VERSION 0.5.1
  GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/libfdf.git"
  GIT_TAG "master"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/libfdf
  )

