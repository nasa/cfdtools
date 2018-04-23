add_executable(cones_of_sight
  cones_of_sight.f90
)
target_link_libraries(cones_of_sight global)
target_link_libraries(cones_of_sight
  xyzq_io
  interp3d
  numodules
  adt
  linsys
  optlib
  geomlib
  searchlib
  intrinsics
)
install(FILES README DESTINATION doc RENAME cones_of_sight)
install(TARGETS cones_of_sight DESTINATION ${CMAKE_INSTALL_BINDIR})
