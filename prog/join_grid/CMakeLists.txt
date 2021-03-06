add_executable(join_grid
  join_grid.f
)
target_link_libraries(join_grid cfdtools_global)
install(FILES README DESTINATION doc RENAME join_grid)
install(TARGETS join_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
