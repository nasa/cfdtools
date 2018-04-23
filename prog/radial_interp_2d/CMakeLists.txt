add_executable(radial_interp_2d
  radial_interp_2d.f90
)
target_link_libraries(radial_interp_2d global)
target_link_libraries(radial_interp_2d
  xyzq_io
  interplib
  numodules
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME radial_interp_2d)
install(TARGETS radial_interp_2d DESTINATION ${CMAKE_INSTALL_BINDIR})
