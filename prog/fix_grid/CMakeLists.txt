add_executable(fix_grid
  fix_grid.f90
)
target_link_libraries(fix_grid
  cfdtools_global
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_interp3d
  cfdtools_searchlib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_xyzq_io
)
install(TARGETS fix_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
