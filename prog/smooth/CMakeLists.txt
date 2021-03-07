add_executable(smooth
  dummies.f
  second.f
  smooth.f
)
target_link_libraries(smooth cfdtools_global)
target_link_libraries(smooth
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_interplib
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_obsolete
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_forsythe
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME smooth)
install(TARGETS smooth DESTINATION ${CMAKE_INSTALL_BINDIR})
