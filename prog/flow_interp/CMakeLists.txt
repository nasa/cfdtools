add_executable(flow_interp
  flow_interp.f90
)
target_link_libraries(flow_interp cfdtools_global)
target_link_libraries(flow_interp
  cfdtools_xyzq_io
  cfdtools_adt
  cfdtools_linsys
  cfdtools_kdtree
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME flow_interp)
install(TARGETS flow_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
