add_executable(heat_shield
  second.f
  heat_shield.f90
  capsules.f90
)
target_link_libraries(heat_shield
  global
  integrate1d
  interplib
  gridlib
  geomlib
  numodules
  linsys
  searchlib
  npoptlib
  blaslib
  blas2lib
  intrinsics
)
install(FILES README DESTINATION doc RENAME heat_shield)
install(TARGETS heat_shield DESTINATION ${CMAKE_INSTALL_BINDIR})
