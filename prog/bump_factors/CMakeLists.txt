add_executable(bump_factors
  bump_factors.f90
)
target_link_libraries(bump_factors global)
target_link_libraries(bump_factors
  numodules
  progtools
  adt
  linsys
  tecplot_io
  geomlib
)
install(FILES README DESTINATION doc RENAME bump_factors)
install(TARGETS bump_factors DESTINATION ${CMAKE_INSTALL_BINDIR})
