add_executable(combine
  combine.f
)
target_link_libraries(combine
  cfdtools_global
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
)
install(TARGETS combine DESTINATION ${CMAKE_INSTALL_BINDIR})
