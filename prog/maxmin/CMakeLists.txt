add_executable(maxmin
  dummies.f
  maxmin.f
  second.f
)
target_link_libraries(maxmin cfdtools_global)
target_link_libraries(maxmin
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
install(FILES README DESTINATION doc RENAME maxmin)
install(TARGETS maxmin DESTINATION ${CMAKE_INSTALL_BINDIR})
