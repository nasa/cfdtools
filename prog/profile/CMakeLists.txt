add_executable(profile
  profile.f
)
target_link_libraries(profile cfdtools_global)
target_link_libraries(profile
  cfdtools_interp3d
  cfdtools_interplib
  cfdtools_numodules
  cfdtools_integratelib
  cfdtools_progtools
  cfdtools_interp2d
  cfdtools_linsys
  cfdtools_obsolete
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_forsythe
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME profile)
install(TARGETS profile DESTINATION ${CMAKE_INSTALL_BINDIR})
