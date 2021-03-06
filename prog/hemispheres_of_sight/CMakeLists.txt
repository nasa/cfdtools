add_executable(hemispheres_of_sight
  hemispheres_of_sight.f90
)
target_link_libraries(hemispheres_of_sight cfdtools_global)
target_link_libraries(hemispheres_of_sight
  triangulation_io
  interp3d
  numodules
  progtools
  adt
  linsys
  surface_patch_utilities
  eigenlib
  optlib
  geomlib
  gridlib
  ugridlib
  xyzq_io
  searchlib
  geomlib
  intrinsics
)
install(FILES README DESTINATION doc RENAME hemispheres_of_sight)
install(TARGETS hemispheres_of_sight DESTINATION ${CMAKE_INSTALL_BINDIR})
