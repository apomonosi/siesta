add_library(${PROJECT_NAME}.aux_contrib_ap OBJECT
       makebox.f  fillbox.f inver3.f  hit.f
       opnout.f  test_ani.f test_md.f test_mdc.f wraxsf1.f wraxsf2.f
       read_xv.f  intpl04.f displa.f read_ev.f itochar.f w_arrow.f
       w_movie.f
)

target_link_libraries(${PROJECT_NAME}.aux_contrib_ap
   PRIVATE
    ${PROJECT_NAME}.libsys
    ${PROJECT_NAME}.libunits)

siesta_add_executable(${PROJECT_NAME}.eig2bxsf
   eig2bxsf.f
   $<TARGET_OBJECTS:${PROJECT_NAME}.aux_contrib_ap>
)
target_link_libraries(${PROJECT_NAME}.eig2bxsf
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

siesta_add_executable(${PROJECT_NAME}.xv2xsf
   xv2xsf.f
   $<TARGET_OBJECTS:${PROJECT_NAME}.aux_contrib_ap>
)
target_link_libraries(${PROJECT_NAME}.xv2xsf
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

siesta_add_executable(${PROJECT_NAME}.md2axsf
   md2axsf.f
   $<TARGET_OBJECTS:${PROJECT_NAME}.aux_contrib_ap>
)
target_link_libraries(${PROJECT_NAME}.md2axsf
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

siesta_add_executable(${PROJECT_NAME}.rho2xsf
   rho2xsf.f
   $<TARGET_OBJECTS:${PROJECT_NAME}.aux_contrib_ap>
)
target_link_libraries(${PROJECT_NAME}.rho2xsf
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

siesta_add_executable(${PROJECT_NAME}.vib2xsf
   vib2xsf.f
   $<TARGET_OBJECTS:${PROJECT_NAME}.aux_contrib_ap>
)
target_link_libraries(${PROJECT_NAME}.vib2xsf
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

siesta_add_executable(${PROJECT_NAME}.fmpdos fmpdos.f)

siesta_add_executable(${PROJECT_NAME}.xv2vesta
   xv2vesta.f
   $<TARGET_OBJECTS:${PROJECT_NAME}.aux_contrib_ap>
)

target_link_libraries(${PROJECT_NAME}.xv2vesta
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

siesta_add_executable(${PROJECT_NAME}.vib2vesta
   vib2vesta.f
   $<TARGET_OBJECTS:${PROJECT_NAME}.aux_contrib_ap>
)
target_link_libraries(${PROJECT_NAME}.vib2vesta
                      PRIVATE
		      ${PROJECT_NAME}.libsys
		      ${PROJECT_NAME}.libunits)

if( SIESTA_INSTALL )
  install(
    TARGETS
      ${PROJECT_NAME}.eig2bxsf
      ${PROJECT_NAME}.xv2xsf
      ${PROJECT_NAME}.md2axsf
      ${PROJECT_NAME}.rho2xsf
      ${PROJECT_NAME}.vib2xsf
      ${PROJECT_NAME}.fmpdos
      ${PROJECT_NAME}.xv2vesta
      ${PROJECT_NAME}.vib2vesta
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

add_subdirectory(Tests)