add_executable(integrate
  integrate.f
)
target_link_libraries(integrate
  global
  geomlib
  integrate1d
  linsys
  progtools
  numodules
  obsolete
  searchlib
  table_io
)
install(TARGETS integrate DESTINATION ${CMAKE_INSTALL_BINDIR})
