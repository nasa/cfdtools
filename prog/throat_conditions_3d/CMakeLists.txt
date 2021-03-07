add_executable(throat_conditions_3d
  throat_conditions_3d.f90
)
target_link_libraries(throat_conditions_3d
  cfdtools_global
  cfdtools_adt
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_interplib
  cfdtools_interp2d
  cfdtools_interp3d
  cfdtools_integratelib
  cfdtools_intrinsics
  cfdtools_lewis
  cfdtools_linsys
  cfdtools_numodules
  cfdtools_optlib
  cfdtools_progtools
  cfdtools_searchlib
  cfdtools_xyzq_io
)
install(FILES README DESTINATION doc RENAME throat_conditions_3d)
install(TARGETS throat_conditions_3d DESTINATION ${CMAKE_INSTALL_BINDIR})
