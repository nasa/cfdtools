add_executable(surface_interp
  surface_interp.f90
)
target_link_libraries(surface_interp cfdtools_global)
target_link_libraries(surface_interp
  progtools
  linsys
  tecplot_io
  tecplot_adt
  geomlib
)
install(FILES README DESTINATION doc RENAME surface_interp)
install(TARGETS surface_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
