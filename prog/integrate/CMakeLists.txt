add_executable(integrate
  integrate.f
)
target_link_libraries(integrate
  cfdtools_global
  geomlib
  integratelib
  linsys
  progtools
  numodules
  obsolete
  searchlib
  table_io
)
install(TARGETS integrate DESTINATION ${CMAKE_INSTALL_BINDIR})
