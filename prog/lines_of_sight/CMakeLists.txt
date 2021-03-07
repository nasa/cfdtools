add_executable(lines_of_sight
  lines_of_sight.f90
)
target_link_libraries(lines_of_sight cfdtools_global)
target_link_libraries(lines_of_sight
  cfdtools_xyzq_io
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_adt
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME lines_of_sight)
install(TARGETS lines_of_sight DESTINATION ${CMAKE_INSTALL_BINDIR})
