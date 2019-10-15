add_executable(integrate
  integrate.f
)
target_link_libraries(integrate
  global
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
