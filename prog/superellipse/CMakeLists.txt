add_executable(superellipse
  superellipse.f
)
target_link_libraries(superellipse
  cfdtools_global
  cfdtools_gridlib
  cfdtools_intrinsics
  cfdtools_integratelib
  cfdtools_interplib
  cfdtools_progtools
  cfdtools_numodules
  cfdtools_searchlib
)
install(TARGETS superellipse DESTINATION ${CMAKE_INSTALL_BINDIR})
