add_executable(fix_filler
  fix_filler.f90
)
target_link_libraries(fix_filler
  cfdtools_global
  cfdtools_adt
  cfdtools_xyzq_io
  cfdtools_geomlib
  cfdtools_gridlib
  cfdtools_optlib
  cfdtools_interp3d
  cfdtools_searchlib
  cfdtools_numodules
)
install(TARGETS fix_filler DESTINATION ${CMAKE_INSTALL_BINDIR})
