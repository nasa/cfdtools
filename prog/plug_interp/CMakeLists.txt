add_executable(plug_interp
  plug_interp.f90
)
target_link_libraries(plug_interp
  cfdtools_global
  cfdtools_geomlib
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
  cfdtools_xyzq_io
)
install(TARGETS plug_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
