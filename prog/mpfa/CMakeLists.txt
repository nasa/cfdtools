add_executable(mpfa
  MEDLI_analysis.f90
)
target_link_libraries(mpfa
  global
  interplib
  searchlib
  numodules
)
install(TARGETS mpfa DESTINATION ${CMAKE_INSTALL_BINDIR})
