add_executable(bsprofile
  bsprofile.f
  getdis.f
)
target_link_libraries(bsprofile global)
target_link_libraries(bsprofile
  interplib
  integrate
  numodules
  progtools
  dt_nurbs
  optlib
  gridlib
  geomlib
  searchlib
  aa_nurbs
)
install(FILES README DESTINATION doc RENAME bsprofile)
install(TARGETS bsprofile DESTINATION ${CMAKE_INSTALL_BINDIR})
