add_executable(cones_of_sight
  cones_of_sight.f90
)
target_link_libraries(cones_of_sight cfdtools_global)
target_link_libraries(cones_of_sight
  cfdtools_xyzq_io
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_adt
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_geomlib
  cfdtools_searchlib
  cfdtools_intrinsics
)
install(FILES README DESTINATION doc RENAME cones_of_sight)
install(TARGETS cones_of_sight DESTINATION ${CMAKE_INSTALL_BINDIR})
