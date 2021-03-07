add_executable(regrid1
  regrid1.f90
)
target_link_libraries(regrid1
  cfdtools_global
  cfdtools_gridlib
  cfdtools_interp2d
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_linsys
  cfdtools_geomlib
  cfdtools_numodules
)
install(TARGETS regrid1 DESTINATION ${CMAKE_INSTALL_BINDIR})
