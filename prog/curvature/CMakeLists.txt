add_executable(curvature
  curvature.f
)
target_link_libraries(curvature cfdtools_global)
target_link_libraries(curvature
  interp3d
  interplib
  numodules
  progtools
  linsys
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME curvature)
install(TARGETS curvature DESTINATION ${CMAKE_INSTALL_BINDIR})
