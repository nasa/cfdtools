add_executable(surface_pad
  surface_pad.f90
)
target_link_libraries(surface_pad cfdtools_global)
target_link_libraries(surface_pad
  xyzq_io
  interp3d
  numodules
  progtools
  surface_patch_utilities
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME surface_pad)
install(TARGETS surface_pad DESTINATION ${CMAKE_INSTALL_BINDIR})
