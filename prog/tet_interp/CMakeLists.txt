add_executable(tet_interp
  tet_interp.f90
)
target_link_libraries(tet_interp cfdtools_global)
target_link_libraries(tet_interp
  cfdtools_triangulation_io
  cfdtools_interp3d
  cfdtools_progtools
  cfdtools_adt
  cfdtools_linsys
  cfdtools_eigenlib
  cfdtools_gridlib
  cfdtools_xyzq_io
  cfdtools_geomlib
)
install(TARGETS tet_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
