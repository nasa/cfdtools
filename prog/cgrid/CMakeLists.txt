add_executable(cgrid
  cgrid.f
)
target_link_libraries(cgrid
  global
  gridlib
  interplib
  searchlib
  optlib
  geomlib
  numodules
)
install(TARGETS cgrid DESTINATION ${CMAKE_INSTALL_BINDIR})
