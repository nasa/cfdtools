add_executable(hb_grid
  hb_grid.f
  smooth1d.f
  capsules.f90
)
target_link_libraries(hb_grid
  cfdtools_global
  cfdtools_gridlib
  cfdtools_interplib
  cfdtools_interp2d
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_optlib
  cfdtools_linsys
  cfdtools_searchlib
  cfdtools_numodules
  cfdtools_intrinsics
)
install(FILES README DESTINATION doc RENAME hb_grid)
install(TARGETS hb_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
