add_executable(regularize
  regularize.f
)
target_link_libraries(regularize
  cfdtools_global
  cfdtools_progtools
  cfdtools_interplib
  cfdtools_searchlib
  cfdtools_numodules
)
install(TARGETS regularize DESTINATION ${CMAKE_INSTALL_BINDIR})
