add_executable(line_grid
  line_grid.f90
)
target_link_libraries(line_grid global)
install(FILES README DESTINATION doc RENAME line_grid)
install(TARGETS line_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
