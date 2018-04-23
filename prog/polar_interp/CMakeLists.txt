add_executable(polar_interp
  polar_interp.f90
)
target_link_libraries(polar_interp global)
target_link_libraries(polar_interp
  numodules
  interplib
  gridlib
  geomlib
  searchlib
)
install(FILES README DESTINATION doc RENAME polar_interp)
install(TARGETS polar_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
