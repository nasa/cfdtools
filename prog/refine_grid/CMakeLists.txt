add_executable(refine_grid
  densify_grid_block.f90
  refine_grid.f90
)
target_link_libraries(refine_grid global)
target_link_libraries(refine_grid
  interp3d
  interplib
  numodules
  linsys
  optlib
  gridlib
  xyzq_io
  searchlib
  geomlib
)
install(FILES README DESTINATION doc RENAME refine_grid)
install(TARGETS refine_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
