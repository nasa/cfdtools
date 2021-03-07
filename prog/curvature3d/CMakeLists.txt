add_executable(curvature3d
  curvature3d.f
)
target_link_libraries(curvature3d cfdtools_global)
target_link_libraries(curvature3d
  cfdtools_interp3d
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME curvature3d)
install(TARGETS curvature3d DESTINATION ${CMAKE_INSTALL_BINDIR})
