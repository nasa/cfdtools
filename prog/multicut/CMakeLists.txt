add_executable(multicut
  multicut.f
)
target_link_libraries(multicut
  cfdtools_global
  cfdtools_numodules
  cfdtools_ugridlib
)
install(TARGETS multicut DESTINATION ${CMAKE_INSTALL_BINDIR})
