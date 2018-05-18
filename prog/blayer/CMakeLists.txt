add_executable(blayer
  blayer.f90
)
target_link_libraries(blayer global)
target_link_libraries(blayer
  interplib
  integrate1d
  numodules
  progtools
  tecplot_io
  tecplot_xyzq_io
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME blayer)
install(TARGETS blayer DESTINATION ${CMAKE_INSTALL_BINDIR})
