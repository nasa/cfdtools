add_executable(radial_interp_2d
  radial_interp_2d.f90
)
target_link_libraries(radial_interp_2d cfdtools_global)
target_link_libraries(radial_interp_2d
  cfdtools_xyzq_io
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME radial_interp_2d)
install(TARGETS radial_interp_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
