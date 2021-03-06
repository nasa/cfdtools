add_executable(interp_1d
  interp_1d.f
)
target_link_libraries(interp_1d cfdtools_global)
target_link_libraries(interp_1d
  numodules
  searchlib
  linsys
  interplib
)
install(FILES README DESTINATION doc RENAME interp_1d)
install(TARGETS interp_1d DESTINATION ${CMAKE_INSTALL_BINDIR})
