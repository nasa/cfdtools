add_executable(fix_filler
  fix_filler.f90
)
target_link_libraries(fix_filler
  cfdtools_global
  adt
  xyzq_io
  geomlib
  gridlib
  optlib
  interp3d
  searchlib
  numodules
)
install(TARGETS fix_filler DESTINATION ${CMAKE_INSTALL_BINDIR})
