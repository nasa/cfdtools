add_executable(tri_to_quad
  tri_to_quad.f90
)
target_link_libraries(tri_to_quad cfdtools_global)
target_link_libraries(tri_to_quad
  cfdtools_triangulation_io
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_tecplot_io
  cfdtools_tecplot_adt
  cfdtools_tecplot_xyzq_io
  cfdtools_ugridlib
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME tri_to_quad)
install(TARGETS tri_to_quad DESTINATION ${CMAKE_INSTALL_BINDIR})
