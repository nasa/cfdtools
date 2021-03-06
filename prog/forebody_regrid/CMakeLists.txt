add_executable(forebody_regrid
  forebody_regrid.f90
  grid_block_structure.f90
)
target_link_libraries(forebody_regrid cfdtools_global)
target_link_libraries(forebody_regrid
  xyzq_io
  interp3d
  numodules
  progtools
  adt
  linsys
  surface_patch_utilities
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME forebody_regrid)
install(TARGETS forebody_regrid DESTINATION ${CMAKE_INSTALL_BINDIR})
