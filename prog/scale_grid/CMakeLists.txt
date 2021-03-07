add_executable(scale_grid
  scale_grid.f
  second.f
)
target_link_libraries(scale_grid cfdtools_global cfdtools_progtools cfdtools_geomlib cfdtools_xyzq_io)
install(TARGETS scale_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
