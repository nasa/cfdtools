add_executable(lines_of_sight
  lines_of_sight.f90
)
target_link_libraries(lines_of_sight cfdtools_global)
target_link_libraries(lines_of_sight
  xyzq_io
  interp3d
  numodules
  progtools
  adt
  linsys
  optlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME lines_of_sight)
install(TARGETS lines_of_sight DESTINATION ${CMAKE_INSTALL_BINDIR})
