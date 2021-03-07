add_executable(surface_diffs
  surface_diffs.f90
)
target_link_libraries(surface_diffs cfdtools_global)
target_link_libraries(surface_diffs
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_tecplot_io
  cfdtools_tecplot_adt
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME surface_diffs)
install(TARGETS surface_diffs DESTINATION ${CMAKE_INSTALL_BINDIR})
