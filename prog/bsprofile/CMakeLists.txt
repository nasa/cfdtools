add_executable(bsprofile
  bsprofile.f
  getdis.f
)
target_link_libraries(bsprofile cfdtools_global)
target_link_libraries(bsprofile
  cfdtools_interplib
  cfdtools_integratelib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_dt_nurbs
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
  cfdtools_aa_nurbs
)
install(FILES README DESTINATION doc RENAME bsprofile)
install(TARGETS bsprofile DESTINATION ${CMAKE_INSTALL_BINDIR})
