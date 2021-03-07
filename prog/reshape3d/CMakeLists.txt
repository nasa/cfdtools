add_executable(reshape3d
  reshape3d.f90
)
target_link_libraries(reshape3d cfdtools_global)
target_link_libraries(reshape3d
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_gridlib
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME reshape3d)
install(TARGETS reshape3d DESTINATION ${CMAKE_INSTALL_BINDIR})
