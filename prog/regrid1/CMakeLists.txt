add_executable(regrid1
  regrid1.f90
)
target_link_libraries(regrid1
  global
  gridlib
  interp2d
  interplib
  searchlib
  linsys
  geomlib
  numodules
)
install(TARGETS regrid1 DESTINATION ${CMAKE_INSTALL_BINDIR})
