add_executable(surface_interp
  surface_interp.f90
)
target_link_libraries(surface_interp cfdtools_global)
target_link_libraries(surface_interp
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_tecplot_io
  cfdtools_tecplot_adt
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME surface_interp)
install(TARGETS surface_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
