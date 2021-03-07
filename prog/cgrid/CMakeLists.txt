add_executable(cgrid
  cgrid.f
)
target_link_libraries(cgrid
  cfdtools_global
  cfdtools_gridlib
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_optlib
  cfdtools_geomlib
  cfdtools_numodules
)
install(TARGETS cgrid DESTINATION ${CMAKE_INSTALL_BINDIR})
