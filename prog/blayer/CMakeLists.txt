add_executable(blayer
  blayer.f90
)
target_link_libraries(blayer cfdtools_global)
target_link_libraries(blayer
  cfdtools_interplib
  cfdtools_integratelib
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_tecplot_io
  cfdtools_tecplot_xyzq_io
  cfdtools_linsys
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME blayer)
install(TARGETS blayer DESTINATION ${CMAKE_INSTALL_BINDIR})
