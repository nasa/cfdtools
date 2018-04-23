add_executable(split_grid
  split_grid.f
)
target_link_libraries(split_grid global)
target_link_libraries(split_grid
  xyzq_io
)
install(FILES README DESTINATION doc RENAME split_grid)
install(TARGETS split_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
