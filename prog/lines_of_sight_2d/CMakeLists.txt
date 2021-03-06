add_executable(lines_of_sight_2d
  lines_of_sight_2d.f90
)
target_link_libraries(lines_of_sight_2d cfdtools_global)
target_link_libraries(lines_of_sight_2d
  xyzq_io
  interplib
  numodules
  progtools
  adt
  intrinsics
  interp2d
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME lines_of_sight_2d)
install(TARGETS lines_of_sight_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
