add_executable(surface_cuts
  surface_cuts.f
)
target_link_libraries(surface_cuts
  cfdtools_global
  gridlib
  geomlib
  interplib
  integratelib
  intrinsics
  numodules
  searchlib
  triangulation_io
  ugridlib
)
install(TARGETS surface_cuts DESTINATION ${CMAKE_INSTALL_BINDIR})
