add_executable(curvature3d
  curvature3d.f
)
target_link_libraries(curvature3d global)
target_link_libraries(curvature3d
  interp3d
  interplib
  numodules
  progtools
  linsys
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME curvature3d)
install(TARGETS curvature3d DESTINATION ${CMAKE_INSTALL_BINDIR})
