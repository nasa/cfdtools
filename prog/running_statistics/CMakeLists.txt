add_executable(running_statistics
  running_statistics.f90
)
target_link_libraries(running_statistics global numodules)
install(TARGETS running_statistics DESTINATION ${CMAKE_INSTALL_BINDIR})
