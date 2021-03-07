add_executable(scan_grid
  scan_grid.f90
)
target_link_libraries(scan_grid cfdtools_global)
target_link_libraries(scan_grid
  cfdtools_gridlib
  cfdtools_progtools
  cfdtools_geomlib
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME scan_grid)
install(TARGETS scan_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
