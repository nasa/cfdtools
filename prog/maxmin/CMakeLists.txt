add_executable(maxmin
  dummies.f
  maxmin.f
  second.f
)
target_link_libraries(maxmin global)
target_link_libraries(maxmin
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
install(FILES README DESTINATION doc RENAME maxmin)
install(TARGETS maxmin DESTINATION ${CMAKE_INSTALL_BINDIR})
