add_executable(thin_grid
  thin_grid.f90
)
target_link_libraries(thin_grid global)
target_link_libraries(thin_grid
  xyzq_io
)
install(FILES README DESTINATION doc RENAME thin_grid)
install(TARGETS thin_grid DESTINATION ${CMAKE_INSTALL_BINDIR})
