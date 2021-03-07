add_executable(adjust_grid
  adjust_grid.f90
)
target_link_libraries(adjust_grid cfdtools_global)
target_link_libraries(adjust_grid
  cfdtools_grid_block_structure
  cfdtools_gridlib
  cfdtools_progtools
  cfdtools_geomlib
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME adjust_grid)
install(TARGETS adjust_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
