add_executable(p3d2tec
  p3d2tec.f90
)
target_link_libraries(p3d2tec cfdtools_global)
target_link_libraries(p3d2tec
  cfdtools_progtools
  cfdtools_tecplot_io
  cfdtools_tecplot_xyzq_io
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME p3d2tec)
install(TARGETS p3d2tec DESTINATION ${CMAKE_INSTALL_BINDIR})
