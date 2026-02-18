#
# There are several targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file". Common files are grouped in a pseudo-library object.
#
# This is a major annoyance.
#
set(top_srcdir "${PROJECT_SOURCE_DIR}/Src")

add_library(${PROJECT_NAME}.coop_common_objs OBJECT
        ${top_srcdir}/alloc.F90
        ${top_srcdir}/parallel.F
        ${top_srcdir}/m_io.f
        ${top_srcdir}/moreParallelSubs.F90
        ${top_srcdir}/m_getopts.f90
        ${top_srcdir}/precision.F
        io.f
        subs.f90
        main_vars.f90
        io_hs.f90
        orbital_set.f90
        read_curves.f90
        local_timer.f90
)
target_link_libraries(${PROJECT_NAME}.coop_common_objs
  PUBLIC
    ${PROJECT_NAME}.libunits
    ${PROJECT_NAME}.mpi_macros
  PRIVATE
    ${PROJECT_NAME}.libsys
  )

set(sources_fat
        fat.f90
)
set(sources_mprop
        mprop.f90
)
set(sources_spin_texture
        spin_texture.f90
)

set(sources_read_spin_texture
        read_spin_texture.f90
)

set(sources_dm_creator
        iodm_netcdf.F90
        write_dm.f
        dm_creator.F90
)

siesta_add_executable(${PROJECT_NAME}.fat
   ${sources_fat}
)
target_link_libraries(${PROJECT_NAME}.fat
  PRIVATE
   ${PROJECT_NAME}.coop_common_objs
   ${PROJECT_NAME}.libsys
   )

siesta_add_executable(${PROJECT_NAME}.mprop
   ${sources_mprop}
)
target_link_libraries(${PROJECT_NAME}.mprop
   PRIVATE
   ${PROJECT_NAME}.coop_common_objs
   ${PROJECT_NAME}.libsys
   )

siesta_add_executable(${PROJECT_NAME}.spin_texture
   ${sources_spin_texture}
)
target_link_libraries(${PROJECT_NAME}.spin_texture
   PRIVATE
   ${PROJECT_NAME}.coop_common_objs
   ${PROJECT_NAME}.libsys
   )

siesta_add_executable(${PROJECT_NAME}.read_spin_texture
   ${sources_read_spin_texture}
)

siesta_add_executable(${PROJECT_NAME}.dm_creator
   ${sources_dm_creator}
)

target_link_libraries(${PROJECT_NAME}.dm_creator
   PRIVATE
   ${PROJECT_NAME}.coop_common_objs
   ${PROJECT_NAME}.libsys
   )

if(SIESTA_WITH_NETCDF)
  target_link_libraries(${PROJECT_NAME}.dm_creator
   PRIVATE
   NetCDF::NetCDF_Fortran
 )
  target_compile_definitions(${PROJECT_NAME}.dm_creator
    PRIVATE
    CDF
  )
endif()



if( SIESTA_INSTALL )
  install(
    TARGETS
      ${PROJECT_NAME}.fat
      ${PROJECT_NAME}.mprop
      ${PROJECT_NAME}.spin_texture
      ${PROJECT_NAME}.dm_creator
      ${PROJECT_NAME}.read_spin_texture
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
