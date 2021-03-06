add_executable(fix_grid
  fix_grid.f90
)
target_link_libraries(fix_grid
  cfdtools_global
  geomlib
  gridlib
  interp3d
  searchlib
  numodules
  progtools
  xyzq_io
)
install(TARGETS fix_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
