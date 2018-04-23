add_executable(profile
  profile.f
)
target_link_libraries(profile global)
target_link_libraries(profile
  interp3d
  interplib
  numodules
  integrate
  progtools
  interp2d
  linsys
  obsolete
  optlib
  gridlib
  forsythe
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME profile)
install(TARGETS profile DESTINATION ${CMAKE_INSTALL_BINDIR})
