add_executable(multicut
  multicut.f
)
target_link_libraries(multicut
  global
  numodules
  ugridlib
)
install(TARGETS multicut DESTINATION ${CMAKE_INSTALL_BINDIR})
