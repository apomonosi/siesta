#
message(STATUS "Proceeding to declare variables for on-the-fly libgridxc compilation...")

foreach( opt
         MPI GRID_SP LIBXC
       )
       
  set(WITH_${opt} ${SIESTA_WITH_${opt}})
  message(STATUS ".... Libgridxc: WITH_${opt}=${WITH_${opt}}")

  # This is the new name in (forthcoming) versions of LibGridXC
  set(LIBGRIDXC_WITH_${opt} ${SIESTA_WITH_${opt}})

endforeach()

