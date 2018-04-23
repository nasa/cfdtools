add_executable(tet_interp
  tet_interp.f90
)
target_link_libraries(tet_interp global)
target_link_libraries(tet_interp
  triangulation_io
  interp3d
  progtools
  adt
  linsys
  eigenlib
  gridlib
  xyzq_io
  geomlib
)
install(TARGETS tet_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
