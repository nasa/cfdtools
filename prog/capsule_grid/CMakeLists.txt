add_executable(capsule_grid
  capsule_grid.f90
  capsules.f90
  grid_block_structure.f90
)
target_link_libraries(capsule_grid global)
target_link_libraries(capsule_grid
  interp3d
  interplib
  numodules
  progtools
  adt
  linsys
  surface_patch_utilities
  gridlib
  optlib
  xyzq_io
  searchlib
  geomlib
  intrinsics
)
install(FILES README DESTINATION doc RENAME capsule_grid)
install(TARGETS capsule_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
