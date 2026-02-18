  # Set w90 variable to compile subproject

  message(STATUS "Proceeding to declare variables for on-the-fly wrapper-wannier90 compilation...")

  set(WITH_MPI ${SIESTA_WITH_MPI})
  message(STATUS ".... wrapper-wannier90: WITH_MPI=${WITH_MPI}")

