add_executable(compress2d
  compress2d.f90
)
target_link_libraries(compress2d cfdtools_global)
target_link_libraries(compress2d
  cfdtools_xyzq_io
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
  cfdtools_intrinsics
  cfdtools_interp3d
  cfdtools_interp2d
)
install(FILES README DESTINATION doc RENAME compress2d)
install(TARGETS compress2d DESTINATION ${CMAKE_INSTALL_BINDIR})
