add_executable(distribute
  distribute.f
  shockgrid.f
)
target_link_libraries(distribute cfdtools_global)
target_link_libraries(distribute
  cfdtools_numodules
  cfdtools_interplib
  cfdtools_progtools
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME distribute)
install(TARGETS distribute DESTINATION ${CMAKE_INSTALL_BINDIR})
