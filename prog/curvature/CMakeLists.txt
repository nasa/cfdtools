add_executable(curvature
  curvature.f
)
target_link_libraries(curvature cfdtools_global)
target_link_libraries(curvature
  cfdtools_interp3d
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME curvature)
install(TARGETS curvature DESTINATION ${CMAKE_INSTALL_BINDIR})
