include(SiestaFindPackage)

Siesta_find_package(toml-f
  REQUIRED
  GIT_REPOSITORY "https://github.com/toml-f/toml-f.git"
  GIT_TAG "v0.2.4"
  SOURCE_DIR "${DFTD3_SOURCE_ROOT}/toml-f"
  )
