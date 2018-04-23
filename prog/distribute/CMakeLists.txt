add_executable(distribute
  distribute.f
  shockgrid.f
)
target_link_libraries(distribute global)
target_link_libraries(distribute
  numodules
  interplib
  progtools
  optlib
  gridlib
  searchlib
)
install(FILES README DESTINATION doc RENAME distribute)
install(TARGETS distribute DESTINATION ${CMAKE_INSTALL_BINDIR})
