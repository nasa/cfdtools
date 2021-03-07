add_executable(cross_sections
  cross_sections.f90
)
target_link_libraries(cross_sections cfdtools_global)
target_link_libraries(cross_sections
  cfdtools_interplib
  cfdtools_integratelib
  cfdtools_intrinsics
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_tecplot_io
  cfdtools_tecplot_xyzq_io
  cfdtools_ugridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME cross_sections)
install(TARGETS cross_sections DESTINATION ${CMAKE_INSTALL_BINDIR})
