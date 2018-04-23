add_executable(radial_interp
  radial_interp.f90
)
target_link_libraries(radial_interp global)
target_link_libraries(radial_interp
  xyzq_io
  interp3d
  numodules
  adt
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME radial_interp)
install(TARGETS radial_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
