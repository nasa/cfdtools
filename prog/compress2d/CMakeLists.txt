add_executable(compress2d
  compress2d.f90
)
target_link_libraries(compress2d cfdtools_global)
target_link_libraries(compress2d
  xyzq_io
  interplib
  numodules
  progtools
  linsys
  optlib
  gridlib
  geomlib
  searchlib
  intrinsics
  interp3d
  interp2d
)
install(FILES README DESTINATION doc RENAME compress2d)
install(TARGETS compress2d DESTINATION ${CMAKE_INSTALL_BINDIR})
