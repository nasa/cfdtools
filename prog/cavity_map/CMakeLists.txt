add_executable(cavity_map
  cavity_map.f90
  rotation_matrix.f90
)
target_link_libraries(cavity_map global)
target_link_libraries(cavity_map
  interp3d
  interplib
  numodules
  progtools
  adt
  tecplot_io
  linsys
  eigenlib
  surface_patch_utilities
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME cavity_map)
install(TARGETS cavity_map DESTINATION ${CMAKE_INSTALL_BINDIR})
