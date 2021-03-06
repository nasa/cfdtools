add_executable(aerosurf
  aerosurf.f
  second.f
  surfer.f
  synch.f
)
target_link_libraries(aerosurf cfdtools_global)
target_link_libraries(aerosurf
  interp3d
  interplib
  numodules
  progtools
  interp2d
  linsys
  optlib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME aerosurf)
install(TARGETS aerosurf DESTINATION ${CMAKE_INSTALL_BINDIR})
