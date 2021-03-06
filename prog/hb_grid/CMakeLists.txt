add_executable(hb_grid
  hb_grid.f
  smooth1d.f
  capsules.f90
)
target_link_libraries(hb_grid
  cfdtools_global
  gridlib
  interplib
  interp2d
  geomlib
  gridlib
  optlib
  linsys
  searchlib
  numodules
  intrinsics
)
install(FILES README DESTINATION doc RENAME hb_grid)
install(TARGETS hb_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
