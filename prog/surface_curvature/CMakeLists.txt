add_executable(surface_curvature
  surface_curvature.f90
)
target_link_libraries(surface_curvature cfdtools_global)
target_link_libraries(surface_curvature
  numodules
  progtools
  tecplot_io
  tecplot_gridlib
  tecplot_xyzq_io
  geomlib
)
install(FILES README DESTINATION doc RENAME surface_curvature)
install(TARGETS surface_curvature DESTINATION ${CMAKE_INSTALL_BINDIR})
