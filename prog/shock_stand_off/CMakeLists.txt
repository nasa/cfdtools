add_executable(shock_stand_off
  shock_stand_off.f90
)
target_link_libraries(shock_stand_off cfdtools_global)
target_link_libraries(shock_stand_off
  cfdtools_xyzq_io
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_adt
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME shock_stand_off)
install(TARGETS shock_stand_off DESTINATION ${CMAKE_INSTALL_BINDIR})
