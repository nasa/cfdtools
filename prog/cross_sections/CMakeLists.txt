add_executable(cross_sections
  cross_sections.f90
)
target_link_libraries(cross_sections global)
target_link_libraries(cross_sections
  interplib
  integrate1d
  intrinsics
  numodules
  progtools
  tecplot_io
  tecplot_xyzq_io
  ugridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME cross_sections)
install(TARGETS cross_sections DESTINATION ${CMAKE_INSTALL_BINDIR})
