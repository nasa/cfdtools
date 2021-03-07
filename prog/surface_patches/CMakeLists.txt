add_executable(surface_patches
  surface_patches.f90
)
target_link_libraries(surface_patches cfdtools_global)
target_link_libraries(surface_patches
  cfdtools_xyzq_io
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_surface_patch_utilities
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME surface_patches)
install(TARGETS surface_patches DESTINATION ${CMAKE_INSTALL_BINDIR})
