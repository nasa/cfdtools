add_executable(flow_interp_2d
  flow_interp_2d.f90
)
target_link_libraries(flow_interp_2d cfdtools_global)
target_link_libraries(flow_interp_2d
  cfdtools_xyzq_io
  cfdtools_progtools
  cfdtools_adt
  cfdtools_linsys
  cfdtools_gridlib
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME flow_interp_2d)
install(TARGETS flow_interp_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
