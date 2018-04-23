add_executable(flow_interp
  flow_interp.f90
)
target_link_libraries(flow_interp global)
target_link_libraries(flow_interp
  xyzq_io
  adt
  linsys
  kdtree
  geomlib
)
install(FILES README DESTINATION doc RENAME flow_interp)
install(TARGETS flow_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
