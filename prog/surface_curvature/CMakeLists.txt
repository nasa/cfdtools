add_executable(surface_curvature
  surface_curvature.f90
)
target_link_libraries(surface_curvature cfdtools_global)
target_link_libraries(surface_curvature
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_tecplot_io
  cfdtools_tecplot_gridlib
  cfdtools_tecplot_xyzq_io
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME surface_curvature)
install(TARGETS surface_curvature DESTINATION ${CMAKE_INSTALL_BINDIR})
