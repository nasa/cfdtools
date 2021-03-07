add_executable(check_boundary
  check_boundary.f90
)
target_link_libraries(check_boundary cfdtools_global cfdtools_xyzq_io)
install(TARGETS check_boundary DESTINATION ${CMAKE_INSTALL_BINDIR})
