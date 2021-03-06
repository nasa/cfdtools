add_executable(cgrid
  cgrid.f
)
target_link_libraries(cgrid
  cfdtools_global
  gridlib
  interplib
  searchlib
  optlib
  geomlib
  numodules
)
install(TARGETS cgrid DESTINATION ${CMAKE_INSTALL_BINDIR})
