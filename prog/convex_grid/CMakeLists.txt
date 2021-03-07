add_executable(convex_grid
  convex_grid.f90
)
target_link_libraries(convex_grid cfdtools_global)
target_link_libraries(convex_grid
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_adt
  cfdtools_linsys
  cfdtools_surface_patch_utilities
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_xyzq_io
  cfdtools_searchlib
  cfdtools_geomlib
)
install(TARGETS convex_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
