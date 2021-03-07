add_executable(aerosurf
  aerosurf.f
  second.f
  surfer.f
  synch.f
)
target_link_libraries(aerosurf cfdtools_global)
target_link_libraries(aerosurf
  cfdtools_interp3d
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_interp2d
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME aerosurf)
install(TARGETS aerosurf DESTINATION ${CMAKE_INSTALL_BINDIR})
