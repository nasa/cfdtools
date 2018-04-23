add_executable(scan_grid
  scan_grid.f90
)
target_link_libraries(scan_grid global)
target_link_libraries(scan_grid
  gridlib
  progtools
  geomlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME scan_grid)
install(TARGETS scan_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
