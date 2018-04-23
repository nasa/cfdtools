add_executable(smooth2d
  dummies.f
  rdtable.f
  second.f
  smooth2d.f
)
target_link_libraries(smooth2d global)
target_link_libraries(smooth2d
  interplib
  numodules
  progtools
  interp2d
  linsys
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME smooth2d)
install(TARGETS smooth2d DESTINATION ${CMAKE_INSTALL_BINDIR})
