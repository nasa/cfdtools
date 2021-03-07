add_executable(tri_to_tri
  tri_to_tri.f90
)
target_link_libraries(tri_to_tri cfdtools_global)
target_link_libraries(tri_to_tri
  cfdtools_triangulation_io
  cfdtools_xyzq_io
  cfdtools_progtools
  cfdtools_adt
  cfdtools_linsys
  cfdtools_ugridlib
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME tri_to_tri)
install(TARGETS tri_to_tri DESTINATION ${CMAKE_INSTALL_BINDIR})
