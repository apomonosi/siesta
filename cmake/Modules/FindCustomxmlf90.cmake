include(SiestaFindPackage)

Siesta_find_package(xmlf90
  REQUIRED
  GIT_REPOSITORY "https://gitlab.com/siesta-project/libraries/xmlf90.git"
  GIT_TAG "master"
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/External/xmlf90
  )

