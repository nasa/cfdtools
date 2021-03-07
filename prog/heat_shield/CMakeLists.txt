add_executable(heat_shield
  second.f
  heat_shield.f90
  capsules.f90
)
target_link_libraries(heat_shield
  cfdtools_global
  cfdtools_integratelib
  cfdtools_interplib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_numodules
  cfdtools_linsys
  cfdtools_searchlib
  cfdtools_npoptlib
  cfdtools_blaslib
  cfdtools_blas2lib
  cfdtools_intrinsics
)
install(FILES README DESTINATION doc RENAME heat_shield)
install(TARGETS heat_shield DESTINATION ${CMAKE_INSTALL_BINDIR})
