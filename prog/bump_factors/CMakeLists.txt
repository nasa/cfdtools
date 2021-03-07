add_executable(bump_factors
  bump_factors.f90
)
target_link_libraries(bump_factors cfdtools_global)
target_link_libraries(bump_factors
  cfdtools_geomlib
  cfdtools_linsys
  cfdtools_progtools
  cfdtools_numodules
  cfdtools_tecplot_adt
  cfdtools_tecplot_io
)
install(FILES README DESTINATION doc RENAME bump_factors)
install(TARGETS bump_factors DESTINATION ${CMAKE_INSTALL_BINDIR})
