add_executable(reshape3d
  reshape3d.f90
)
target_link_libraries(reshape3d global)
target_link_libraries(reshape3d
  numodules
  progtools
  gridlib
  geomlib
)
install(FILES README DESTINATION doc RENAME reshape3d)
install(TARGETS reshape3d DESTINATION ${CMAKE_INSTALL_BINDIR})
