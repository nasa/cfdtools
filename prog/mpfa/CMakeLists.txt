add_executable(mpfa
  MEDLI_analysis.f90
)
target_link_libraries(mpfa
  cfdtools_global
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
)
install(TARGETS mpfa DESTINATION ${CMAKE_INSTALL_BINDIR})
