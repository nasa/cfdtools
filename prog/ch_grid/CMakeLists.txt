add_executable(ch_grid
  ch_grid.f
  second.f
)
target_link_libraries(ch_grid global)
target_link_libraries(ch_grid
  interp3d
  interplib
  integrate
  intrinsics
  numodules
  progtools
  interp2d
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME ch_grid)
install(TARGETS ch_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
