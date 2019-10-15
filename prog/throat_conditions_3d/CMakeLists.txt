add_executable(throat_conditions_3d
  throat_conditions_3d.f90
)
target_link_libraries(throat_conditions_3d
  global
  adt
  geomlib
  gridlib
  interplib
  interp2d
  interp3d
  integratelib
  intrinsics
  lewis
  linsys
  numodules
  optlib
  progtools
  searchlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME throat_conditions_3d)
install(TARGETS throat_conditions_3d DESTINATION ${CMAKE_INSTALL_BINDIR})
