add_executable(regularize
  regularize.f
)
target_link_libraries(regularize
  cfdtools_global
  progtools
  interplib
  searchlib
  numodules
)
install(TARGETS regularize DESTINATION ${CMAKE_INSTALL_BINDIR})
