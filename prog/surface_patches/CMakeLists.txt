add_executable(surface_patches
  surface_patches.f90
)
target_link_libraries(surface_patches cfdtools_global)
target_link_libraries(surface_patches
  xyzq_io
  interplib
  numodules
  progtools
  surface_patch_utilities
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME surface_patches)
install(TARGETS surface_patches DESTINATION ${CMAKE_INSTALL_BINDIR})
