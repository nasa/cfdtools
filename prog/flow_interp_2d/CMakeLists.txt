add_executable(flow_interp_2d
  flow_interp_2d.f90
)
target_link_libraries(flow_interp_2d cfdtools_global)
target_link_libraries(flow_interp_2d
  xyzq_io
  progtools
  adt
  linsys
  gridlib
  geomlib
)
install(FILES README DESTINATION doc RENAME flow_interp_2d)
install(TARGETS flow_interp_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
