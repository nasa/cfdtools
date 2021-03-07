add_executable(surface_pad
  surface_pad.f90
)
target_link_libraries(surface_pad cfdtools_global)
target_link_libraries(surface_pad
  cfdtools_xyzq_io
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_surface_patch_utilities
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME surface_pad)
install(TARGETS surface_pad DESTINATION ${CMAKE_INSTALL_BINDIR})
