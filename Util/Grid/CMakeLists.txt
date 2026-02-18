#
# There are several targets in this directory. The 'OBJECT library' idiom
# was needed to avoid ninja complaining about "two paths for the same
# .mod file"
#
set(top_src_dir "${PROJECT_SOURCE_DIR}/Src")

# To avoid "multiple rules" in Ninja
# Maybe turn into global-utility targets

siesta_add_library(${PROJECT_NAME}.aux_opts
  NO_LINKER_FLAGS
  OBJECT

  ${top_src_dir}/m_getopts.f90
)

siesta_add_library(${PROJECT_NAME}.aux_gridfunc
  NO_LINKER_FLAGS
  OBJECT

  m_gridfunc.F90
)

#-----------------------------------------
siesta_add_executable(${PROJECT_NAME}.g2c_ng
   ${top_src_dir}/reclat.f
   m_struct.f90
   local_die.f90
   g2c_ng.f
)

# For future enabling of NetCDF in m_gridfunc
target_link_libraries(${PROJECT_NAME}.g2c_ng
  PRIVATE
    $<$<BOOL:${SIESTA_WITH_NETCDF}>:NetCDF::NetCDF_Fortran>
    ${PROJECT_NAME}.aux_opts
    ${PROJECT_NAME}.aux_gridfunc
    ${PROJECT_NAME}.libunits
  )
target_compile_definitions(${PROJECT_NAME}.g2c_ng
  PRIVATE $<$<BOOL:${SIESTA_WITH_NETCDF}>:CDF> )
		     
#---------------------------------------------------------------------

siesta_add_executable(${PROJECT_NAME}.v_info v_info.f90)

# For future enabling of NetCDF in m_gridfunc
target_link_libraries(${PROJECT_NAME}.v_info
  PRIVATE
    $<$<BOOL:${SIESTA_WITH_NETCDF}>:NetCDF::NetCDF_Fortran>
    ${PROJECT_NAME}.aux_gridfunc
    ${PROJECT_NAME}.libunits
  )
target_compile_definitions(${PROJECT_NAME}.v_info  PRIVATE $<$<BOOL:${SIESTA_WITH_NETCDF}>:CDF> )
		     

#---------------------------------------------------------------------
siesta_add_executable(${PROJECT_NAME}.grid2cube grid2cube.f)
siesta_add_executable(${PROJECT_NAME}.grid2val grid2val.F)
siesta_add_executable(${PROJECT_NAME}.grid_rotate grid_rotate.F90)
siesta_add_executable(${PROJECT_NAME}.grid_supercell grid_supercell.F90)

if( SIESTA_INSTALL )
  install(
    TARGETS
      ${PROJECT_NAME}.g2c_ng
      ${PROJECT_NAME}.v_info
      ${PROJECT_NAME}.grid2cube
      ${PROJECT_NAME}.grid2val
      ${PROJECT_NAME}.grid_rotate
      ${PROJECT_NAME}.grid_supercell
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif()
#---------------------------------------------------------------------

# Run executables with NetCDF dependencies
if (SIESTA_WITH_NETCDF)

siesta_add_executable(${PROJECT_NAME}.cdf2grid  cdf2grid.F90)

target_link_libraries(${PROJECT_NAME}.cdf2grid PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(${PROJECT_NAME}.cdf2grid  PRIVATE CDF )
		     
#---------------------------------------------------------------------
siesta_add_executable(${PROJECT_NAME}.grid2cdf  grid2cdf.F90)

target_link_libraries(${PROJECT_NAME}.grid2cdf PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(${PROJECT_NAME}.grid2cdf  PRIVATE CDF )
		     
#-----------------------------------------
siesta_add_executable(${PROJECT_NAME}.cdf2xsf  cdf2xsf.F90)

target_link_libraries(${PROJECT_NAME}.cdf2xsf PRIVATE NetCDF::NetCDF_Fortran )
target_compile_definitions(${PROJECT_NAME}.cdf2xsf  PRIVATE CDF )
		     
#-----------------------------------------

# To avoid "multiple rules" in Ninja
# (could add OBJECT qualifier to avoid the lib appearing)
siesta_add_library(${PROJECT_NAME}.aux_m_grid  m_grid.F90)
target_compile_definitions(${PROJECT_NAME}.aux_m_grid  PRIVATE CDF )
target_link_libraries(${PROJECT_NAME}.aux_m_grid PUBLIC NetCDF::NetCDF_Fortran )


siesta_add_executable(${PROJECT_NAME}.cdf_diff cdf_diff.F90)

target_link_libraries(${PROJECT_NAME}.cdf_diff
  PRIVATE
    ${PROJECT_NAME}.aux_m_grid
    NetCDF::NetCDF_Fortran
)
target_compile_definitions(${PROJECT_NAME}.cdf_diff  PRIVATE CDF )
		     
#-----------------------------------------
siesta_add_executable(${PROJECT_NAME}.cdf_get_cell cdf_get_cell.F90)

target_link_libraries(${PROJECT_NAME}.cdf_get_cell
  PRIVATE
    ${PROJECT_NAME}.aux_m_grid
    NetCDF::NetCDF_Fortran
)
target_compile_definitions(${PROJECT_NAME}.cdf_get_cell  PRIVATE CDF )
		     
#-----------------------------------------

# To avoid "multiple rules" in Ninja
siesta_add_library(${PROJECT_NAME}.aux_fft
  NO_LINKER_FLAGS
  OBJECT
    ${top_src_dir}/m_fft_gpfa.F
)
target_link_libraries(${PROJECT_NAME}.aux_fft PRIVATE ${PROJECT_NAME}.libsys)

siesta_add_executable(${PROJECT_NAME}.cdf_laplacian
   ${top_src_dir}/reclat.f
   fft3d.f
   local_die.f90
   cdf_laplacian.F90
)

target_link_libraries(${PROJECT_NAME}.cdf_laplacian
  PRIVATE
    ${PROJECT_NAME}.aux_fft
    ${PROJECT_NAME}.aux_opts
    NetCDF::NetCDF_Fortran
)
target_compile_definitions(${PROJECT_NAME}.cdf_laplacian  PRIVATE CDF )
		     
#-----------------------------------------

siesta_add_executable(${PROJECT_NAME}.cdf_fft
   ${top_src_dir}/reclat.f
   fft3d.f
   local_die.f90
   cdf_fft.F90
)

target_link_libraries(${PROJECT_NAME}.cdf_fft
  PRIVATE
    ${PROJECT_NAME}.aux_fft
    NetCDF::NetCDF_Fortran
)
target_compile_definitions(${PROJECT_NAME}.cdf_fft  PRIVATE CDF )
		     
if( SIESTA_INSTALL )
  install(TARGETS
    ${PROJECT_NAME}.cdf2grid
    ${PROJECT_NAME}.grid2cdf
    ${PROJECT_NAME}.cdf2xsf
    ${PROJECT_NAME}.cdf_diff
    ${PROJECT_NAME}.cdf_get_cell
    ${PROJECT_NAME}.cdf_laplacian
    ${PROJECT_NAME}.cdf_fft
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} )
endif()

endif(SIESTA_WITH_NETCDF)
