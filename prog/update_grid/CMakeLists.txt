add_executable(update_grid
  update_grid.f90
)
target_link_libraries(update_grid cfdtools_global)
target_link_libraries(update_grid
  cfdtools_progtools
  cfdtools_geomlib
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME update_grid)
install(TARGETS update_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
