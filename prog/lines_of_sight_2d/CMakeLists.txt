add_executable(lines_of_sight_2d
  lines_of_sight_2d.f90
)
target_link_libraries(lines_of_sight_2d cfdtools_global)
target_link_libraries(lines_of_sight_2d
  cfdtools_xyzq_io
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_adt
  cfdtools_intrinsics
  cfdtools_interp2d
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME lines_of_sight_2d)
install(TARGETS lines_of_sight_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
