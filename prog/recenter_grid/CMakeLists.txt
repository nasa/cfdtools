add_executable(recenter_grid
  recenter_grid.f90
)
target_link_libraries(recenter_grid cfdtools_global xyzq_io)
install(TARGETS recenter_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
