add_executable(blayer
  blayer.f90
  extract_blayer_data.f90
)
target_link_libraries(blayer global)
target_link_libraries(blayer
  xyzq_io
  interplib
  integrate
  numodules
  progtools
  tecplot_io
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME blayer)
install(TARGETS blayer DESTINATION ${CMAKE_INSTALL_BINDIR})
