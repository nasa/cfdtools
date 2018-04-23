add_executable(gen1d
  gen1d.f90
)
target_link_libraries(gen1d global)
target_link_libraries(gen1d
  numodules
  interplib
  progtools
  linsys
  special_functions
  optlib
  gridlib
  geomlib
)
install(FILES README DESTINATION doc RENAME gen1d)
install(TARGETS gen1d DESTINATION ${CMAKE_INSTALL_BINDIR})
