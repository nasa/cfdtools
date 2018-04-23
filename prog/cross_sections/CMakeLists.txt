add_executable(cross_sections
  cross_sections.f90
)
target_link_libraries(cross_sections global)
target_link_libraries(cross_sections
  xyzq_io
  interplib
  integrate
  numodules
  progtools
  tecplot_io
  ugridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME cross_sections)
install(TARGETS cross_sections DESTINATION ${CMAKE_INSTALL_BINDIR})
