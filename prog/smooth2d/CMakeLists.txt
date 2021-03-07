add_executable(smooth2d
  dummies.f
  rdtable.f
  second.f
  smooth2d.f
)
target_link_libraries(smooth2d cfdtools_global)
target_link_libraries(smooth2d
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_interp2d
  cfdtools_linsys
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME smooth2d)
install(TARGETS smooth2d DESTINATION ${CMAKE_INSTALL_BINDIR})
