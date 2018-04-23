add_executable(tri_to_tri
  tri_to_tri.f90
)
target_link_libraries(tri_to_tri global)
target_link_libraries(tri_to_tri
  triangulation_io
  xyzq_io
  progtools
  adt
  linsys
  ugridlib
  geomlib
)
install(FILES README DESTINATION doc RENAME tri_to_tri)
install(TARGETS tri_to_tri DESTINATION ${CMAKE_INSTALL_BINDIR})
