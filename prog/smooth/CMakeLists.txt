add_executable(smooth
  dummies.f
  second.f
  smooth.f
)
target_link_libraries(smooth global)
target_link_libraries(smooth
  interp3d
  numodules
  interplib
  progtools
  linsys
  obsolete
  optlib
  gridlib
  forsythe
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME smooth)
install(TARGETS smooth DESTINATION ${CMAKE_INSTALL_BINDIR})
