add_executable(tri_to_quad
  tri_to_quad.f90
)
target_link_libraries(tri_to_quad global)
target_link_libraries(tri_to_quad
  triangulation_io
  progtools
  linsys
  tecplot_io
  tecplot_adt
  tecplot_xyzq_io
  ugridlib
  geomlib
)
install(FILES README DESTINATION doc RENAME tri_to_quad)
install(TARGETS tri_to_quad DESTINATION ${CMAKE_INSTALL_BINDIR})
