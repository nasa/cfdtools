add_executable(shock_stand_off
  shock_stand_off.f90
)
target_link_libraries(shock_stand_off global)
target_link_libraries(shock_stand_off
  xyzq_io
  interp3d
  numodules
  adt
  linsys
  optlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME shock_stand_off)
install(TARGETS shock_stand_off DESTINATION ${CMAKE_INSTALL_BINDIR})
