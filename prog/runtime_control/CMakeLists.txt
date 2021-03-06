add_executable(runtime_control
  runtime.f90
)
target_link_libraries(runtime_control cfdtools_global progtools)
install(TARGETS runtime_control DESTINATION ${CMAKE_INSTALL_BINDIR})
