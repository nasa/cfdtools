add_executable(gen1d
  gen1d.f90
)
target_link_libraries(gen1d cfdtools_global)
target_link_libraries(gen1d
  cfdtools_numodules
  cfdtools_interplib
  cfdtools_progtools
  cfdtools_linsys
  cfdtools_special_functions
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
)
install(FILES README DESTINATION doc RENAME gen1d)
install(TARGETS gen1d DESTINATION ${CMAKE_INSTALL_BINDIR})
