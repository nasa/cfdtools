add_executable(convex_grid
  convex_grid.f90
)
target_link_libraries(convex_grid cfdtools_global)
target_link_libraries(convex_grid
  interp3d
  numodules
  progtools
  adt
  linsys
  surface_patch_utilities
  optlib
  gridlib
  xyzq_io
  searchlib
  geomlib
)
install(TARGETS convex_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
