add_executable(p3d2tec
  p3d2tec.f90
)
target_link_libraries(p3d2tec global)
target_link_libraries(p3d2tec
  progtools
  tecplot_io
  geomlib
  xyzq_io
)
install(FILES README DESTINATION doc RENAME p3d2tec)
install(TARGETS p3d2tec DESTINATION ${CMAKE_INSTALL_BINDIR})
