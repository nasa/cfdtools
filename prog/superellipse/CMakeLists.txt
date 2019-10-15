add_executable(superellipse
  superellipse.f
)
target_link_libraries(superellipse
  global
  gridlib
  intrinsics
  integratelib
  interplib
  progtools
  numodules
  searchlib
)
install(TARGETS superellipse DESTINATION ${CMAKE_INSTALL_BINDIR})
