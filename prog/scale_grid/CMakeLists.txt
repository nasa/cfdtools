add_executable(scale_grid
  scale_grid.f
  second.f
)
target_link_libraries(scale_grid global progtools geomlib xyzq_io)
install(TARGETS scale_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
