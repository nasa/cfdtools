add_executable(tri_to_quad
  Tec_dummy.f90
  Tecplot_io.f90
  tri_to_quad.f90
)
target_link_libraries(tri_to_quad global)
target_link_libraries(tri_to_quad
  triangulation_io
  xyzq_io
  progtools
  adt
  linsys
  ugridlib
  geomlib
)
install(FILES README DESTINATION doc RENAME tri_to_quad)
install(TARGETS tri_to_quad DESTINATION ${CMAKE_INSTALL_BINDIR})
