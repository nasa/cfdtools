add_executable(plug_interp
  plug_interp.f90
)
target_link_libraries(plug_interp
  cfdtools_global
  geomlib
  interplib
  searchlib
  numodules
  xyzq_io
)
install(TARGETS plug_interp DESTINATION ${CMAKE_INSTALL_BINDIR})
