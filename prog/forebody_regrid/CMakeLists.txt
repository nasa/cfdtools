add_executable(forebody_regrid
  forebody_regrid.f90
  grid_block_structure.f90
)
target_link_libraries(forebody_regrid cfdtools_global)
target_link_libraries(forebody_regrid
  cfdtools_xyzq_io
  cfdtools_interp3d
  cfdtools_numodules
  cfdtools_progtools
  cfdtools_adt
  cfdtools_linsys
  cfdtools_surface_patch_utilities
  cfdtools_optlib
  cfdtools_gridlib
  cfdtools_geomlib
  cfdtools_searchlib
)
install(FILES README DESTINATION doc RENAME forebody_regrid)
install(TARGETS forebody_regrid DESTINATION ${CMAKE_INSTALL_BINDIR})
