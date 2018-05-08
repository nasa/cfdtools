add_executable(bump_factors
  bump_factors.f90
)
target_link_libraries(bump_factors global)
target_link_libraries(bump_factors
  geomlib
  linsys
  progtools
  numodules
  tecplot_adt
  tecplot_io
)
install(FILES README DESTINATION doc RENAME bump_factors)
install(TARGETS bump_factors DESTINATION ${CMAKE_INSTALL_BINDIR})
