add_executable(surface_interp_2d
  surface_interp_2d.f90
)
target_link_libraries(surface_interp_2d cfdtools_global)
target_link_libraries(surface_interp_2d
  progtools
  linsys
  tecplot_io
  tecplot_adt
  geomlib
)
install(FILES README DESTINATION doc RENAME surface_interp_2d)
install(TARGETS surface_interp_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
