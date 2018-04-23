add_executable(surface_diffs
  surface_diffs.f90
)
target_link_libraries(surface_diffs global)
target_link_libraries(surface_diffs
  progtools
  adt
  linsys
  tecplot_io
  geomlib
)
install(FILES README DESTINATION doc RENAME surface_diffs)
install(TARGETS surface_diffs DESTINATION ${CMAKE_INSTALL_BINDIR})
